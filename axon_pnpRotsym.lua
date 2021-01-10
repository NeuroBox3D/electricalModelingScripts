--------------------------------------------------------------------------------
-- 3D PNP on axon with Ranvier nodes.                                         --
-- The implementation relies on perfect rotational symmetry to reduce the     --
-- discretization to a 2d problem.                                            --
-- Parallel execution is recommended.                                         --
--                                                                            --
-- Author: Markus Breit                                                       --
-- Date:   2018-09-03                                                         --
--------------------------------------------------------------------------------

-- load pre-implemented lua functions
ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")

AssertPluginsLoaded({"neuro_collection", "nernst_planck", "ConvectionDiffusion", "Parmetis"}) 

-- speed up evaluation of lua functions by c program
EnableLUA2C(true)

-- choose algebra
InitUG(2, AlgebraType("CPU", 5))


------------------------------------
-- process command line arguments --
------------------------------------
-- grid
grid = util.GetParam("-grid", "../grids/axon5_2d_pnp.ugx")
nRanvier = util.GetParam("-nRanv", 5)

-- time stepping
numTimeSteps =  util.GetParamNumber("-nSteps", 500, "number of timesteps")
dt = util.GetParamNumber("-dt", 2e-5, "time step in seconds")  -- in s

-- spatial discretization
numRefs = util.GetParamNumber("-numRefs", 2, "number of refinements")
numAnisoRefs = util.GetParamNumber("-numAnisoRefs", 3, "number of anisotropic refinements")
numMemRefs = util.GetParamNumber("-numMemRefs", 0, "number of (hanging) membrane refinements")

-- choose outfile directory
outDir = util.GetParam("-outName", "test/neuron")
outDir = outDir.."/"

-- generate vtk output? (every modulo-th step)
generateVTKoutput = util.HasParamOption("-vtk")
modulo = util.GetParamNumber("-modulo", 1)

-- specify "-verbose" to output linear solver convergence
verbose	= util.HasParamOption("-verbose")


---------------------------------------
-- problem constants and stimulation --
---------------------------------------
Phi				=	"Phi"
ion_species		=	{   "K"   ,  "Cl"  ,  "Na"  ,   "A" }
valency			=	{    1    ,   -1   ,   1    ,   -1  }
conc_out      	= 	{     4   ,   123  ,  145   ,   26  }   -- in mol/m^3 == mM
conc_in			= 	{   155   ,     5  ,   12   ,  162  }   -- in mol/m^3 == mM
diff_const		=	{  1.96e-9, 2.03e-9, 1.33e-9, 2.0e-9}   -- in m^2/s

F = 96485  -- in C / mol
R = 8.31451  -- in J / (mol K)
T = 298.15  -- in K
vac_perm = 8.8541878e-12  -- in C / (Vm)

rel_perm_mem = 4.0
rel_perm_rest = 80.0
perm_mem = rel_perm_mem * vac_perm
perm_rest = rel_perm_rest * vac_perm

pot_out = 0.0
pot_in = -0.070

-- channel conductances
g_Na = 1.2e+03  -- in C / (V s m^2)
g_K = 3.6e+02  -- in C / (V s m^2)
g_L = 3.0e+00  -- in C / (V s m^2)

-- equilibrium potentials
E_K = R*T/F*math.log(conc_out[1]/conc_in[1])  -- -0.0940  -- in V
E_Na = R*T/F*math.log(conc_out[3]/conc_in[3])  -- 0.0640  -- in V
leak_K = g_K * math.pow(n_inf(1e3*pot_in),4) * (pot_in - E_K) -- 0.03092
leak_Na = g_Na * math.pow(m_inf(1e3*pot_in),3)*h_inf(1e3*pot_in) * (pot_in - E_Na) -- -0.00293
E_L_K = pot_in + 2*leak_K/g_L
E_L_Na = pot_in + 2*leak_Na/g_L

-- additional scaling of equations
pot_fac = 1.0/F


-- dimless parameters --
conc_ref = 100.0
pot_ref = R*T/F
diff_ref = 2e-9
len_ref = 1e-6
time_ref = len_ref*len_ref/diff_ref
flux_ref = diff_ref * conc_ref / len_ref

dimless_endTime = numTimeSteps*dt / time_ref
dimless_dt = dt / time_ref

dimless_conc_out = {}
dimless_conc_in = {}
dimless_diff_const = {}

for i = 1,4 do
	dimless_conc_out[i] = conc_out[i] / conc_ref
	dimless_conc_in[i] = conc_in[i] / conc_ref
	dimless_diff_const[i] = diff_const[i] / diff_ref
end

dimless_pot_out = pot_out / pot_ref
dimless_pot_in = pot_in / pot_ref

dimless_rel_perm_mem = rel_perm_mem * pot_ref / (len_ref*len_ref*conc_ref) * pot_fac
dimless_rel_perm_rest = rel_perm_rest * pot_ref / (len_ref*len_ref*conc_ref) * pot_fac
dimless_perm_mem = dimless_rel_perm_mem * vac_perm
dimless_perm_rest = dimless_rel_perm_rest * vac_perm

dim_scales = {pot_ref, conc_ref, conc_ref, conc_ref, conc_ref}

-- injection current density
injCurrentDensity = 3.0  -- in C / (s m^2)
injectionStart = 0  -- in s
injectionEnd = 1  -- in s

function electrodeCurrent_in(z, r, t, si)
	if t*time_ref >= injectionStart and t*time_ref <= injectionEnd then
		return r*injCurrentDensity / (-F * flux_ref)
	end
	return 0.0
end


-------------------------------
-- approximation space setup --
-------------------------------
neededSubsets = {"in", "out", "myelin", "myelin_in", "myelin_out", "bnd_in", "bnd_out",
				 "ranvier1", "ranvier1_in", "ranvier1_out"}
for i = 2,nRanvier do
	neededSubsets[10+3*(i-2)+1] = "ranvier" .. i
	neededSubsets[10+3*(i-2)+2] = "ranvier" .. i .. "_in"
	neededSubsets[10+3*(i-2)+3] = "ranvier" .. i .. "_out"
end
dom = util.CreateDomain(grid, 0, neededSubsets)

-- rescale domain
if len_ref ~= 1e-6 then  -- domain is saved in units of um
	local sel = Selector(dom:grid())
	SelectDomainElements(sel, true, true, false, false, false)
	ScaleDomain(dom, sel, MakeVec(0, 0, 0), MakeVec(1e-6/len_ref, 1e-6/len_ref, 1.0)) 
end

balancer.partitioner = "parmetis"
balancer.staticProcHierarchy = true
balancer.firstDistLvl = 0
balancer.redistSteps = 3
balancer.firstDistProcs = 48
balancer.redistProcs = 8
balancer.parallelElementThreshold = 24
balancer.ParseParameters()

loadBalancer = balancer.CreateLoadBalancer(dom)
if loadBalancer ~= nil then
	if balancer.partitioner == "parmetis" then
		mu = ManifoldUnificator(dom)
		mu:add_protectable_subsets("myelin_in")
		mu:add_protectable_subsets("myelin_out")
		for i = 1,nRanvier do
			mu:add_protectable_subsets("ranvier"..i.."_in")
			mu:add_protectable_subsets("ranvier"..i.."_out")
		end
		cdgm = ClusteredDualGraphManager()
		cdgm:add_unificator(SiblingUnificator())
		cdgm:add_unificator(mu)
		au = AnisotropyUnificator(dom)
		au:set_threshold_ratio(0.1)
		cdgm:add_unificator(au)
		balancer.defaultPartitioner:set_dual_graph_manager(cdgm)
	end
	balancer.Rebalance(dom, loadBalancer)
	loadBalancer:estimate_distribution_quality()
end

-- name some useful subset collections
ranvier = "ranvier1"
ranvier_in = "ranvier1_in"
ranvier_out = "ranvier1_out"
for i = 2,nRanvier do
	ranvier = ranvier .. ", ranvier" .. i
	ranvier_in = ranvier_in .. ", ranvier" .. i .. "_in"
	ranvier_out = ranvier_out .. ", ranvier" .. i .. "_out"
end
memAll = "myelin, " .. ranvier
memAll_in = "myelin_in, " .. ranvier_in
memAll_out = "myelin_out, " .. ranvier_out
volOut = "out, bnd_out, " .. memAll
volIn = "in, bnd_in, " .. memAll
injection = "ranvier1"

-- create approximation space
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("Phi", "Lagrange", 1)
approxSpace:add_fct(ion_species, "Lagrange", 1)
approxSpace:init_top_surface()

-- anisotropic refinements
refiner = HangingNodeDomainRefiner(dom)
AddShadowCopyAdjuster(refiner)
for i = 1, numAnisoRefs do
	MarkAnisotropicX(refiner, dom, 0.76)
	unmark_ranvier_areas(refiner, approxSpace, ranvier, i <= 3)
	refiner:refine()
	balancer.qualityRecordName = "anisoRef " .. i
	balancer.Rebalance(dom, loadBalancer)
	if loadBalancer ~= nil then
		loadBalancer:estimate_distribution_quality()
	end
end


-- membrane refinements
surfaceSubsets = {"myelin_in", "myelin_out"}
volumeSubsets = {"in", "out"}
for i = 1, nRanvier do
	table.insert(surfaceSubsets, "ranvier"..i.."_in")
	table.insert(surfaceSubsets, "ranvier"..i.."_out")
	table.insert(volumeSubsets, "in")
	table.insert(volumeSubsets, "out")
end

for i = 1, numMemRefs do
	MarkAlongSurface(refiner, dom, surfaceSubsets, volumeSubsets)
	refiner:refine()
	balancer.qualityRecordName = "  memRef " .. i
	balancer.Rebalance(dom, loadBalancer)
	if loadBalancer ~= nil then
		loadBalancer:estimate_distribution_quality()
	end
end


-- isotropic refinements
for i = 1, numRefs do
	MarkGlobal(refiner, dom)
	refiner:refine()
	balancer.qualityRecordName = " globRef " .. i
	balancer.Rebalance(dom, loadBalancer)
	if loadBalancer ~= nil then
		loadBalancer:estimate_distribution_quality()
	end
end
delete(refiner)


-- init and order dof distros AFTER complete refinement and redistribution
approxSpace:init_levels()
if numMemRefs == 0 then
	OrderCuthillMcKee(approxSpace, true)
	--OrderLex(approxSpace, "y")
end

if loadBalancer ~= nil then
	loadBalancer:print_quality_records()
end
print()
print(dom:domain_info():to_string())

--SaveGridHierarchyTransformed(dom:grid(), dom:subset_handler(), outDir.."grid/refined_grid_hierarchy_p" .. ProcRank() .. ".ugx", 0.2)
--SaveParallelGridLayout(dom:grid(), outDir.."grid/parallel_grid_layout_p"..ProcRank()..".ugx", 2.0)


--------------------------
-- discretization setup --
--------------------------
-- create domain disc
domainDisc = DomainDiscretization(approxSpace)

function scaled_perm(z,r,t,si)
	if si == 0 or si == 1 then
		local c = r * dimless_perm_rest
		return  c, 0,
				0, c
	end
	local c = r * dimless_perm_mem
	return  c, 0,
			0, c
end
function scaled_diff_1(z,r,t)
	--LUACompiler:ignore
	local c = r * dimless_diff_const[1]
	return  c, 0,
			0, c
end
function scaled_diff_2(z,r,t)
	--LUACompiler:ignore
	local c = r * dimless_diff_const[2]
	return  c, 0,
			0, c
end
function scaled_diff_3(z,r,t)
	--LUACompiler:ignore
	local c = r * dimless_diff_const[3]
	return  c, 0,
			0, c
end
function scaled_diff_4(z,r,t)
	--LUACompiler:ignore
	local c = r * dimless_diff_const[4]
	return  c, 0,
			0, c
end

function ranvier_flux_K(z, r, t, si)
	--LUACompiler:ignore
	if si >= 10 and si <= 19 then
		local ranvNode = math.floor((si-10)/2)+1
		vm = ranvierData[ranvNode].vm
		n = ranvierData[ranvNode].n
		if si % 2 == 0 then
			return -r*(g_K*n*n*n*n*(vm - E_K) + 0.5*g_L*(vm - E_L_K)) / (valency[1]*F*flux_ref)
		end
		return r*(g_K*n*n*n*n*(vm - E_K) + 0.5*g_L*(vm - E_L_K)) / (valency[1]*F*flux_ref)
	end
	return 0.0
end
function ranvier_flux_Na(z, r, t, si)
	--LUACompiler:ignore
	if si >= 10 and si <= 19 then
		local ranvNode = math.floor((si-10)/2)+1
		vm = ranvierData[ranvNode].vm
		m = ranvierData[ranvNode].m
		h = ranvierData[ranvNode].h
		if si % 2 == 0 then
			return -r*(g_Na*m*m*m*h*(vm - E_Na) + 0.5*g_L*(vm - E_L_Na)) / (valency[3]*F*flux_ref)
		end
		return r*(g_Na*m*m*m*h*(vm - E_Na) + 0.5*g_L*(vm - E_L_Na)) / (valency[3]*F*flux_ref)
	end
	return 0.0
end

function rotsym_scale(z,r,t)
	return r
end
rotsymScale = LuaUserNumber("rotsym_scale")



-- inside and outside: fully coupled Nernst-Planck equations
sourceCharge = ScaleAddLinkerNumber()
potentialDisc = ConvectionDiffusionFV1(Phi, "in, out")
potentialDisc:set_stationary()
potentialDisc:set_diffusion("scaled_perm")
potentialDisc:set_source(sourceCharge)
domainDisc:add(potentialDisc)
	
for i = 1, #ion_species do
	-- compute velocity
	local vel = ScaleAddLinkerVector()
	vel:add(-valency[i]*F/(R*T)*pot_ref*dimless_diff_const[i]*rotsymScale, potentialDisc:gradient())

	-- equation for species
	local speciesDisc = ConvectionDiffusionFV1(ion_species[i], "in, out")
	speciesDisc:set_diffusion("scaled_diff_"..i)
	speciesDisc:set_velocity(vel)
	speciesDisc:set_mass_scale(rotsymScale)
	speciesDisc:set_upwind(FullUpwind())
	domainDisc:add(speciesDisc)
	
	-- contribution to charge source
	sourceCharge:add(valency[i]*F*pot_fac*rotsymScale, speciesDisc:value())			
end

-- membrane: only electro-static problem for potential, ohmic conductivity on Ranvier nodes
estatPotentialDisc = ConvectionDiffusionFV1(Phi, memAll)
estatPotentialDisc:set_stationary()
estatPotentialDisc:set_diffusion("scaled_perm")
domainDisc:add(estatPotentialDisc)

membraneDiscK = UserFluxBoundaryFV1("K", ranvier_in .. ", " .. ranvier_out)
membraneDiscK:set_flux_function("ranvier_flux_K")
domainDisc:add(membraneDiscK)

membraneDiscNa = UserFluxBoundaryFV1("Na", ranvier_in .. ", " .. ranvier_out)
membraneDiscNa:set_flux_function("ranvier_flux_Na")
domainDisc:add(membraneDiscNa)

memSpec = DirichletBoundary()
for i = 1, #ion_species do
	memSpec:add(0.0, ion_species[i], memAll)
end
domainDisc:add(memSpec)

-- upper Dirichlet boundary
diriBnds = DirichletBoundary()
diriBnds:add(dimless_pot_out, Phi, "bnd_out")
for i = 1, #ion_species do
	diriBnds:add(dimless_conc_out[i], ion_species[i], "bnd_out")
end
domainDisc:add(diriBnds)

-- extra Dirichlet bnd (only for stationary case)
statDiriBnds = DirichletBoundary()
statDiriBnds:add(dimless_pot_in, Phi, "bnd_in")
for i = 1, #ion_species do
	statDiriBnds:add(dimless_conc_in[i], ion_species[i], "bnd_in")
end
domainDisc:add(statDiriBnds)

-- constraints for adaptivity
hangingConstraint = SymP1Constraints()--OneSideP1Constraints()
domainDisc:add(hangingConstraint)

-- time discretization
timeDisc = ThetaTimeStep(domainDisc)
timeDisc:set_theta(1.0)


-------------------
-- algebra setup --
-------------------
-- debug writer
dbgWriter = GridFunctionDebugWriter(approxSpace)
dbgWriter:set_vtk_output(false)
dbgWriter:set_base_dir(outDir)

convCheck = ConvCheck()
convCheck:set_maximum_steps(100)
convCheck:set_minimum_defect(1e-30)
convCheck:set_reduction(1e-06)
convCheck:set_verbose(verbose)
--[[
convCheck = CompositeConvCheck(approxSpace, 10, 1e-30, 1e-08)
convCheck:set_component_check(Phi, 1e-30, 1e-08)
convCheck:set_component_check(ion_species, 1e-30, 1e-08)
convCheck:set_verbose(true)
convCheck:set_time_measurement(true)
--]]

-- solver setup --
linearSolver = BiCGStab()
gmg = GeometricMultiGrid(approxSpace)
gmg:set_discretization(timeDisc)
gmg:set_base_level(0)
gmg:set_gathered_base_solver_if_ambiguous(true)
--gmg:set_debug(dbgWriter)

base = SuperLU()
gmg:set_base_solver(base)

smoother = ILU()
if numMemRefs ~= 0 then
	smoother:set_sort(true)  -- neither Cuthill-McKee nor Lex seem to be robust with hanging nodes
end
smoother:set_beta(-0.05)
smoother:enable_consistent_interfaces(true)
gmg:set_smoother(smoother)
gmg:set_smooth_on_surface_rim(true)
gmg:set_cycle_type(1)
gmg:set_num_presmooth(3)
gmg:set_num_postsmooth(3)
gmg:set_rap(true)

linearSolver:set_preconditioner(gmg)
linearSolver:set_convergence_check(convCheck)
--linearSolver:set_debug(dbgWriter)


--- non-linear solver ---
-- convergence check
newtonConvCheck = CompositeConvCheck(approxSpace, 10, 1e-10, 1e-12)
newtonConvCheck:set_component_check(Phi, 1e-08, 1e-12)
newtonConvCheck:set_component_check(ion_species, 1e-09, 1e-12)
newtonConvCheck:set_verbose(true)
newtonConvCheck:set_time_measurement(true)
--newtonConvCheck:set_adaptive(true)

-- Newton solver
newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(linearSolver)
newtonSolver:set_convergence_check(newtonConvCheck)
--newtonSolver:set_debug(dbgWriter)


------------------------------------------------------
--  Initial values  ----------------------------------
------------------------------------------------------

u = GridFunction(approxSpace)
u:set(0.0)
u_scaled = u:clone()

InterpolateInner(dimless_pot_out, u, Phi, "out, bnd_out, meas_out, " .. memAll_out .. ", " .. memAll)
InterpolateInner(dimless_pot_in, u, Phi, memAll_in .. ", in, bnd_in, meas_in")

for i = 1, #ion_species do
	InterpolateInner(dimless_conc_out[i], u, ion_species[i], "out, bnd_out, meas_out, " .. memAll_out .. ", " .. ranvier)
	InterpolateInner(dimless_conc_in[i] , u, ion_species[i], memAll_in .. ", in, bnd_in, meas_in")
end

ranvierData = {}
for r = 1, nRanvier do
	ranvierData[r] = {}
	ranvierData[r].vm = pot_in - pot_out
	ranvierData[r].n = n_inf(1e3*(pot_in - pot_out))
	ranvierData[r].m = m_inf(1e3*(pot_in - pot_out))
	ranvierData[r].h = h_inf(1e3*(pot_in - pot_out))
end


------------------------------------------------------
--  compute start value for time loop  ---------------
------------------------------------------------------
op_stat = AssembledOperator(domainDisc)
op_stat:init()
newtonSolver:init(op_stat)

if newtonSolver:apply(u) == false then
	print("Newton solver failed.")
	exit()
end


-------------------
-- Time stepping --
-------------------
-- initial values
time = 0
step = 0

-- Hodgkin-Huxley mechanism update function
area_ranvier_in = {}
area_ranvier_out = {}
for r = 1, nRanvier do
	area_ranvier_in[r] = Integral(1.0, u, "ranvier"..r.."_in")
	area_ranvier_out[r] = Integral(1.0, u, "ranvier"..r.."_out")
end
function updateRanvier(u, dt)
	for r = 1, nRanvier do
		local vm = 1e3 * pot_ref * (Integral(u, Phi, "ranvier"..r.."_in")/area_ranvier_in[r]
		                            - Integral(u, Phi, "ranvier"..r.."_out")/area_ranvier_out[r])
		local ninf = n_inf(vm)
		local minf = m_inf(vm)
		local hinf = h_inf(vm)
		local taun = tau_n(vm)
		local taum = tau_m(vm)
		local tauh = tau_h(vm)
		
		ranvierData[r].vm = 1e-3*vm
		ranvierData[r].n = ninf + (ranvierData[r].n - ninf) * math.exp(-dt*1e3/taun)
		ranvierData[r].m = minf + (ranvierData[r].m - minf) * math.exp(-dt*1e3/taum)
		ranvierData[r].h = hinf + (ranvierData[r].h - hinf) * math.exp(-dt*1e3/tauh)
	end
end


-- remove stationary Dirichlet bnds on inner bnd (GMG transfer operators must be updated!)
domainDisc:remove(statDiriBnds)
gmg:force_reinit()


-- injection current
elemDiscInjection_in = UserFluxBoundaryFV1(ion_species[4], "ranvier1_in")  -- function order: to, from
elemDiscInjection_in:set_flux_function("electrodeCurrent_in")
domainDisc:add(elemDiscInjection_in)


-- write start solutions
ScaleGF(u_scaled, u, dim_scales)
if generateVTKoutput then
	out = VTKOutput()
	out:print(outDir.."vtk/solution", u_scaled, step, time)
end

-- measure initial potential
take_measurement(u_scaled, time, ranvier_in .. ", " .. ranvier_out, "Phi, K, Na, Cl, A", outDir.."meas/meas")
take_measurement(u_scaled, time, "meas_in, meas_out", "Phi, K, Na, Cl, A", outDir.."meas/bulk")

-- create new grid functions for old value and rhs
uOld = u:clone()

-- store grid function in vector of old solutions
solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, time)

-- init operator
op = AssembledOperator(timeDisc)
newtonSolver:init(op)

-- reset convergence check params
newtonConvCheck:set_component_check(Phi, 1e-13, 1.1)
newtonConvCheck:set_component_check(ion_species, 1e-06*dt, 1e-08)

-- start the time stepping
for step = 1, numTimeSteps do
	print("++++++ TIMESTEP " .. step .. " BEGIN ++++++")
	
	-- prepare step
	updateRanvier(u, dt)
	timeDisc:prepare_step(solTimeSeries, dt/time_ref)

	if not newtonSolver:apply(u) then
		print("Error: Solver for membrane potential problem did not converge.")
		if generateVTKoutput then
			out:write_time_pvd(outDir.."vtk/solution", u)
		end
		exit()
	end
			
	-- update to new time
	time = time + dt
		
	-- plot solution
	ScaleGF(u_scaled, u, dim_scales)
	if generateVTKoutput and step % modulo == 0 then
		out:print(outDir.."vtk/solution", u_scaled, step/modulo, time)
	end
	
	-- measure potential
	take_measurement(u_scaled, time, ranvier_in .. ", " .. ranvier_out, "Phi, K, Na, Cl, A", outDir.."meas/meas")
	take_measurement(u_scaled, time, "meas_in, meas_out", "Phi, K, Na, Cl, A", outDir.."meas/bulk")

	-- update time series
	oldSol = solTimeSeries:oldest()
	VecScaleAssign(oldSol, 1.0, u)
	solTimeSeries:push_discard_oldest(oldSol, time)
	
	print("++++++ TIMESTEP " .. step .. " END ++++++++")
end


-- end timeseries, produce gathering file
if generateVTKoutput then
    out:write_time_pvd(outDir.."vtk/solution", u)
end

