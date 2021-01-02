--------------------------------------------------------------------------------
-- 1D-PNP on cylindrical axon with Ranvier nodes.                             --
--                                                                            --
-- Author: Markus Breit                                                       --
-- Date:   2018-09-04                                                         --
--------------------------------------------------------------------------------

-- load pre-implemented lua functions
ug_load_script("ug_util.lua")

AssertPluginsLoaded({"cable_neuron", "neuro_collection", "nernst_planck", "ConvectionDiffusion"}) 

-- speed up evaluation of lua functions by c program
EnableLUA2C(true)

-- init with dimension and algebra
InitUG(1, AlgebraType("CPU", 5))


------------------------------------
-- process command line arguments --
------------------------------------
-- grid
grid = util.GetParam("-grid", "grids/axon5_1d.ugx")
nRanvier = util.GetParam("-nRanv", 5)
axon_radius = util.GetParamNumber("-rad", 2.5e-7)  -- in m

-- time stepping
numTimeSteps =  util.GetParamNumber("-nSteps", 500, "number of timesteps")
dt = util.GetParamNumber("-dt", 2e-5, "time step in seconds")  -- in s

-- spatial discretization
numRefs = util.GetParamNumber("-numRefs", 4, "number of refinements")

-- choose outfile directory
outDir = util.GetParam("-outName", "test/neuron")
outDir = outDir.."/"

-- generate vtk output? (every modulo-th step)
generateVTKoutput = util.HasParamOption("-vtk")
modulo = util.GetParamNumber("-modulo", 1)

-- specify "-verbose" to output linear solver convergence
verbose	= util.HasParamOption("-verbose")


-----------------------
-- problem constants --
-----------------------
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

--[[
sigma_i = 0.0
sigma_o = 0.0
for i = 1, #ion_species do
	sigma_i = sigma_i + diff_const[i]*valency[i]*F*valency[i]*F/(R*T)*conc_in[i]
	sigma_o = sigma_o + diff_const[i]*valency[i]*F*valency[i]*F/(R*T)*conc_out[i]
end
print("sigma_i = " .. sigma_i)
print("sigma_o = " .. sigma_o)
---]]

-- geometry
radius = axon_radius
mem_thickness_ranvier = radius*math.log(1.0 + 1e-8/radius) --1e-8  -- 10nm
mem_thickness_myelin = radius*math.log(1.0 + 2e-7/radius) --20.0*mem_thickness_ranvier  -- 10 layers of myelin

-- membrane capacitances
total_conc = 0
for i=1,4 do
	total_conc = total_conc + valency[i]*valency[i]*conc_in[i]
end
cap_factor_ranvier = perm_mem / mem_thickness_ranvier / total_conc
cap_factor_myelin = perm_mem / mem_thickness_myelin / total_conc
spec_cap_ranvier = {}  -- in F/m^2
spec_cap_myelin = {}  -- in F/m^2
for i=1,4 do
	spec_cap_ranvier[i] = valency[i]*valency[i]*conc_in[i] * cap_factor_ranvier
	spec_cap_myelin[i] = valency[i]*valency[i]*conc_in[i] * cap_factor_myelin
end

-- channel conductances
g_Na = 1.2e+03  -- in C / (V s m^2)
g_K = 3.6e+02  -- in C / (V s m^2)
g_L = 3.0e+00  -- in C / (V s m^2)

-- equilibrium potentials
E_Na = 0.0640  -- in V
E_K = -0.0940  -- in V
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
cap_ref = conc_ref * len_ref / pot_ref * F
flux_ref = diff_ref * conc_ref / len_ref

dimless_endTime = numTimeSteps*dt / time_ref
dimless_dt = dt / time_ref

dimless_conc_out = {}
dimless_conc_in = {}
dimless_diff_const = {}
dimless_spec_cap_ranvier = {}
dimless_spec_cap_myelin = {}

for i = 1,4 do
	dimless_conc_out[i] = conc_out[i] / conc_ref
	dimless_conc_in[i] = conc_in[i] / conc_ref
	dimless_diff_const[i] = diff_const[i] / diff_ref
	dimless_spec_cap_ranvier[i] = spec_cap_ranvier[i] / cap_ref
	dimless_spec_cap_myelin[i] = spec_cap_myelin[i] / cap_ref
end

dimless_pot_out = pot_out / pot_ref
dimless_pot_in = pot_in / pot_ref

dimless_rel_perm_rest = rel_perm_rest * pot_ref / (len_ref*len_ref*conc_ref) * pot_fac
dimless_perm_rest = dimless_rel_perm_rest * vac_perm

dimless_radius = radius / len_ref

dim_scales = {pot_ref, conc_ref, conc_ref, conc_ref, conc_ref}

print("Reference time        : " .. time_ref)
print("Reference length      : " .. len_ref)
print("Reference potential   : " .. pot_ref)
print("Reference conentration: " .. conc_ref)

-- injection current density
injCurrentDensity = 3.0  -- in C / (s m^2)
injectionStart = 0  -- in s
injectionEnd = 1.0  -- in s

function electrodeCurrent(x, t, si)
	---[[
	-- constant stimulation (windowed)
	if t*time_ref >= injectionStart and t*time_ref <= injectionEnd then
		return 2.0*math.pi*dimless_radius *
			injCurrentDensity / (-F * flux_ref)
	end
	return 0.0
	--]]
	--[[ -- sinusoidal simulation
	return 2.0*math.pi*dimless_radius *
		injCurrentDensity / (-F * flux_ref) *
		0.5*(math.cos(2*math.pi*t*time_ref/0.03)+1)
	--]]
end


-------------------------
-- approximation space --
-------------------------
-- create domain
dom = Domain()
dom:create_additional_subset_handler("projSH")
LoadDomain(dom, grid)

-- check subsets
neededSubsets = {"myelin", "ranvier1"}
for i = 2,nRanvier do
	neededSubsets[4+i] = "ranvier" .. i
end
ug_assert(util.CheckSubsets(dom, neededSubsets), "Something wrong with required subsets.")

-- refine
if numRefs > 0 then
	local refiner = GlobalDomainRefiner(dom)
	for i = 1, numRefs do
		TerminateAbortedRun()
		refiner:refine()
	end
	delete(refiner)
end

-- rescale domain
if len_ref ~= 1 then  -- domain is saved in units of m
	scale_domain(dom, 1.0/len_ref)  -- this scales the positions AND diameter attachment
end

-- name some useful subset collections
ranvier = "ranvier1"
ranvier_vec = {"ranvier1"}
for i = 2,nRanvier do
	ranvier = ranvier .. ", ranvier" .. i
	ranvier_vec[i] = "ranvier" .. i
end
memAll = "myelin, " .. ranvier
injection = "ranvier1"

approxSpace = ApproximationSpace(dom)
approxSpace:add_fct(Phi, "Lagrange", 1)			-- add electric potential
approxSpace:add_fct(ion_species, "Lagrange", 1)	-- add species
approxSpace:init_top_surface()
approxSpace:init_levels()

OrderCuthillMcKee(approxSpace, true)


print(dom:domain_info():to_string())
approxSpace:print_statistic()

--[[
RescaleDomain(dom, len_ref)  -- scale to real size
SaveGridHierarchyTransformed(dom:grid(), dom:subset_handler(), outDir.."refined_grid_hierarchy_p" .. ProcRank() .. ".ugx", 100)
SaveParallelGridLayout(dom:grid(), outDir.."parallel_grid_layout_lv_p"..ProcRank()..".ugx", 100)
RescaleDomain(dom, 1.0/len_ref)  -- scale back to reference length
--]]


--------------------
-- discretization --
--------------------
-- create domain disc
domainDisc = DomainDiscretization(approxSpace)

-- factors by which to multiply all integrals to simulate 3d
dimFacVol = math.pi * dimless_radius * dimless_radius
dimFacSurf = 2.0 * math.pi * dimless_radius


-- potential
sourceCharge = ScaleAddLinkerNumber()
potDisc = ConvectionDiffusionFV1(Phi, memAll)
potDisc:set_stationary()
potDisc:set_diffusion(dimFacVol*dimless_perm_rest)
potDisc:set_source(dimFacVol*sourceCharge)
domainDisc:add(potDisc)

-- species
for i = 1, #ion_species do
	-- drift velocity
	local vel = ScaleAddLinkerVector()
	vel:add(-dimFacVol*dimless_diff_const[i]*valency[i]*F/(R*T)*pot_ref, potDisc:gradient())

	-- capacitive term
	local mass_ranvier = ScaleAddLinkerNumber()
	mass_ranvier:add(dimFacSurf*dimless_spec_cap_ranvier[i]/valency[i], potDisc:value())
	local mass_myelin = ScaleAddLinkerNumber()
	mass_myelin:add(dimFacSurf*dimless_spec_cap_myelin[i]/valency[i], potDisc:value())

	local speciesDisc_ranvier = ConvectionDiffusionFV1(ion_species[i], ranvier)
	speciesDisc_ranvier:set_upwind(FullUpwind())
	speciesDisc_ranvier:set_diffusion(dimFacVol*dimless_diff_const[i])
	speciesDisc_ranvier:set_velocity(vel)
	speciesDisc_ranvier:set_mass(mass_ranvier)
	speciesDisc_ranvier:set_mass_scale(dimFacVol)
	
	local speciesDisc_myelin = ConvectionDiffusionFV1(ion_species[i], "myelin")
	speciesDisc_myelin:set_upwind(FullUpwind())
	speciesDisc_myelin:set_diffusion(dimFacVol*dimless_diff_const[i])
	speciesDisc_myelin:set_velocity(vel)
	speciesDisc_myelin:set_mass(mass_myelin)
	speciesDisc_myelin:set_mass_scale(dimFacVol)
	
	domainDisc:add(speciesDisc_ranvier)
	domainDisc:add(speciesDisc_myelin)
	
	-- contribution to charge source
	sourceCharge:add(valency[i]*F*pot_fac, speciesDisc_myelin:value())			
end

-- HH
hh = HHSpecies({"K", "", "Na", "", "Phi", ""}, ranvier_vec, dom:subset_handler())
hh:set_constant(1, dimless_conc_out[1])
hh:set_constant(3, dimless_conc_out[3])
hh:set_constant(5, 0.0)
hh:set_conductances(g_K, g_Na)
hh:set_reversal_potentials(E_K, E_Na)
hh:set_temperature(T)
hh:use_exact_gating_mode()
hh:set_reference_time(time_ref)
hh:set_scale_inputs({conc_ref, conc_ref, conc_ref, conc_ref, pot_ref, pot_ref})
hh:set_scale_fluxes({1.0/(valency[1]*F*flux_ref), 1.0/(valency[3]*F*flux_ref)})
hhDisc = MembraneTransport1d(ranvier, hh)
hhDisc:set_density_function(1.0)
--hhDisc:set_radius(dimless_radius)  -- read directly from file now
domainDisc:add(hhDisc)

-- leakage
leakageK = OhmicLeakageCharges({ion_species[1], "", Phi, ""})
leakageK:set_constant(1, dimless_conc_out[1])
leakageK:set_constant(3, 0.0)
leakageK:set_conductance(0.5*g_L)
leakageK:set_reversal_potential(E_L_K)
leakageK:set_scale_inputs({conc_ref, conc_ref, pot_ref, pot_ref})
leakageK:set_scale_fluxes({1.0/(valency[1]*F*flux_ref)})
leakageKDisc = MembraneTransport1d(ranvier, leakageK)
leakageKDisc:set_density_function(1.0)
--leakageKDisc:set_radius(dimless_radius)
domainDisc:add(leakageKDisc)

leakageNa = OhmicLeakageCharges({ion_species[3], "", Phi, ""})
leakageNa:set_constant(1, dimless_conc_out[3])
leakageNa:set_constant(3, 0.0)
leakageNa:set_conductance(0.5*g_L)
leakageNa:set_reversal_potential(E_L_Na)
leakageNa:set_scale_inputs({conc_ref, conc_ref, pot_ref, pot_ref})
leakageNa:set_scale_fluxes({1.0/(valency[3]*F*flux_ref)})
leakageNaDisc = MembraneTransport1d(ranvier, leakageNa)
leakageNaDisc:set_density_function(1.0)
--leakageNaDisc:set_radius(dimless_radius)
domainDisc:add(leakageNaDisc)

-- electrode
electrodeDisc = ConvectionDiffusionFV1("A", "ranvier1")
electrodeDisc:set_mass_scale(0.0)
electrodeDisc:set_source("electrodeCurrent")
domainDisc:add(electrodeDisc)


-- time discretization
timeDisc = ThetaTimeStep(domainDisc)
timeDisc:set_theta(1.0) -- 1.0 is implicit Euler

-- create nonlinear operator
op = AssembledOperator(timeDisc)


------------------
-- solver setup --
------------------
convCheck = ConvCheck()
convCheck:set_minimum_defect(1e-50)
convCheck:set_reduction(1e-08)
convCheck:set_verbose(verbose)
convCheck:set_maximum_steps(10000)

linearSolver = LinearSolver()
linearSolver:set_convergence_check(convCheck)
linearSolver:set_preconditioner(ILU())
convCheck:set_maximum_steps(1000)

newtonConvCheck = CompositeConvCheck(approxSpace, 10, 1e-22, 1e-12)
newtonConvCheck:set_component_check(Phi, 1e-05, 1.1)
newtonConvCheck:set_component_check(ion_species, 1e-09*dimless_dt, 1e-08)
newtonConvCheck:set_verbose(true)
newtonConvCheck:set_time_measurement(true)

newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(linearSolver)
newtonSolver:set_convergence_check(newtonConvCheck)


-------------------
-- time stepping --
-------------------
-- initial values
u = GridFunction(approxSpace)
u:set(0.0)

InterpolateInner(dimless_pot_in, u, Phi)
for i = 1, #ion_species do
	InterpolateInner(dimless_conc_in[i], u, ion_species[i])
end

-- init operator
op:init()
newtonSolver:init(op)


dimless_time = 0.0
step = 0

-- write start solutions
u_scaled = u:clone()
ScaleGF(u_scaled, u, dim_scales)
if generateVTKoutput then
	out = VTKOutput()
	out:select_all(true)
	out:print(outDir .."vtk/solution", u_scaled, step, dimless_time*time_ref)
end

-- measure initial potential
take_measurement(u_scaled, 0.0, ranvier, "Phi, K, Na, Cl, A", outDir.."meas/meas")



-- create new grid function for old value
uOld = u:clone()

-- store grid function in vector of old solutions
solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, dimless_time)

while dimless_endTime-dimless_time > 0.001*dimless_dt do
	print("++++++ POINT IN TIME  " .. math.floor((dimless_time+dimless_dt)/dimless_dt+0.5)*dt .. "s  BEGIN ++++++")
	
	-- setup time disc for old solutions and timestep
	timeDisc:prepare_step(solTimeSeries, dimless_dt)
	
	-- apply Newton solver
	if not newtonSolver:apply(u) then 
		print("Error: Solver for membrane potential problem did not converge.")
		if generateVTKoutput then	
			out:write_time_pvd(outDir.."vtk/solution", u)
		end
		exit()
	end
	
	-- update new time
	dimless_time = timeDisc:future_time()
		
	-- plot solution every pstep seconds
	ScaleGF(u_scaled, u, dim_scales)
	if generateVTKoutput and step % modulo == 0 then
		out:print(outDir .."vtk/solution", u_scaled, step/modulo, dimless_time*time_ref)
	end
	
	-- measure potential
	take_measurement(u_scaled, dimless_time*time_ref, ranvier, "Phi, K, Na, Cl, A", outDir.."meas/meas")
	
	-- swap solutions for new time step
	oldestSol = solTimeSeries:oldest()
	VecScaleAssign(oldestSol, 1.0, u)
	solTimeSeries:push_discard_oldest(oldestSol, dimless_time)
	
	step = step + 1
	
	print("++++++ POINT IN TIME  " .. math.floor(dimless_time/dimless_dt+0.5)*dt .. "s  END ++++++++")
end

-- end timeseries, produce gathering file
if generateVTKoutput then
	out:write_time_pvd(outDir .."vtk/solution", u)
end

