--------------------------------------------------------------------------------
-- 3D cable equation on axon with Ranvier nodes.                              --
-- The implementation relies on perfect rotational symmetry to reduce the     --
-- discretization to a 2d problem.                                            --
-- Parallel execution is recommended.                                         --
--                                                                            --
-- Author: Markus Breit                                                       --
-- Date:   2018-09-03                                                         --
--------------------------------------------------------------------------------

ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")

AssertPluginsLoaded({"neuro_collection", "Parmetis"}) 

-- choose algebra
InitUG(2, AlgebraType("CPU", 1))


------------------------------------
-- process command line arguments --
------------------------------------
-- grid
grid = util.GetParam("-grid", "grids/axon5_2d.ugx")
nRanvier = util.GetParam("-nRanv", 5)

-- time stepping
numTimeSteps =  util.GetParamNumber("-nSteps", 2000, "number of timesteps")
dt = util.GetParamNumber("-dt", 0.005, "time step in milliseconds")

-- spatial discretization
numRefs = util.GetParamNumber("-numRefs", 2, "number of refinements")
numAnisoRefs = util.GetParamNumber("-numAnisoRefs", 3, "number of anisotropic refinements")

-- choose outfile directory
outDir = util.GetParam("-outName", "test/neuron")
outDir = outDir.."/"

-- generate vtk output? (every modulo-th step)
generateVTKoutput = util.HasParamOption("-vtk")
modulo = util.GetParamNumber("-modulo", 1) 


---------------------------------------
-- problem constants and stimulation --
---------------------------------------
-- conductivities
sigma_in = 2.4556479647688e3  -- deduced from PNP model, in fC / (ms mV um)
sigma_out = 1.8866027810257e3  -- deduced from PNP model, in fC / (ms mV um)

-- membrane specific capacitance
ranvier_cap = 3.6120456108572e-3  -- 3.54167512e-03  -- in C / (V m^2)
myelin_cap = 2.4101772506797e-4  -- 0.05*ranvier_cap  -- 10 layers with 2 membranes each

-- channel conductances
g_Na = 1.2e-00  -- in fC / (mV ms um^2)
g_K = 3.6e-01  -- in fC / (mV ms um^2)
g_L = 3.0e-03  -- in fC / (mV ms um^2)

-- equilibrium potentials
phi_eq_in = -70.0  -- in mV
phi_eq_out = 0.0  -- in mV
E_Na = 64.0  -- in mV
E_K = -94.0  -- in mV
leak_K = 0.03092
leak_Na = -0.00293
E_L = (phi_eq_in-phi_eq_out) + (leak_K+leak_Na)/g_L

-- injection current density
injCurrentDensity = 3.0  -- in fC / (ms um^2)
injectionStart = 0  -- in ms
injectionEnd = 1000  -- in ms

function electrodeCurrent(z, r, t, si)
	if t >= injectionStart and t <= injectionEnd then
		return 2.0*math.pi*r * injCurrentDensity
	end
	return 0.0
end


-------------------------------
-- approximation space setup --
-------------------------------
neededSubsets = {"in", "out", "myelin", "bnd", "ranvier1"}
for i = 2,nRanvier do
	neededSubsets[4+i] = "ranvier" .. i
end
dom = util.CreateDomain(grid, 0, neededSubsets)

balancer.partitioner = "parmetis"
balancer.staticProcHierarchy = true
balancer.firstDistLvl = 0
balancer.redistSteps = 4
balancer.firstDistProcs = 40
balancer.redistProcs = 64
balancer.parallelElementThreshold = 24
balancer.ParseParameters()

loadBalancer = balancer.CreateLoadBalancer(dom)
if loadBalancer ~= nil then
	if balancer.partitioner == "parmetis" then
		mu = ManifoldUnificator(dom)
		mu:add_protectable_subsets("myelin")
		for i = 1,nRanvier do
			mu:add_protectable_subsets("ranvier"..i)
		end
		au = AnisotropyUnificator(dom)
		au:set_threshold_ratio(0.1)
		cdgm = ClusteredDualGraphManager()
		cdgm:add_unificator(SiblingUnificator())
		cdgm:add_unificator(mu)
		cdgm:add_unificator(au)
		balancer.defaultPartitioner:set_dual_graph_manager(cdgm)
	end
	balancer.Rebalance(dom, loadBalancer)
	loadBalancer:estimate_distribution_quality()
	loadBalancer:print_quality_records()
end

-- name some useful subset collections
ranvier = "ranvier1"
for i = 2,nRanvier do
	ranvier = ranvier .. ", ranvier" .. i
end
memAll = "myelin, " .. ranvier
volOut = "out, bnd, " .. memAll .. ", meas"
volIn = "in, " .. memAll
injection = "ranvier1"

-- create approximation space
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("phiOut", "Lagrange", 1, volOut)
approxSpace:add_fct("phiIn", "Lagrange", 1, volIn)
approxSpace:add_fct("n", "Lagrange", 1, ranvier)
approxSpace:add_fct("m", "Lagrange", 1, ranvier)
approxSpace:add_fct("h", "Lagrange", 1, ranvier)
approxSpace:init_top_surface()
approxSpace:init_levels()

-- anisotropic refinements
refiner = HangingNodeDomainRefiner(dom)
AddShadowCopyAdjuster(refiner)
for i = 1, numAnisoRefs do
	MarkAnisotropicX(refiner, dom, 0.707)
	unmark_ranvier_areas(refiner, approxSpace, ranvier, i <= 3)
	refiner:refine()
	refiner:clear_marks()
	balancer.qualityRecordName = "anisoRef " .. i
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

if loadBalancer ~= nil then
	loadBalancer:print_quality_records()
end
print()
print(dom:domain_info():to_string())

--SaveGridHierarchyTransformed(dom:grid(), dom:subset_handler(), outDir.."grid/refined_grid_hierarchy" .."_p" .. ProcRank() .. ".ugx", 2.0)
--SaveParallelGridLayout(dom:grid(), fileName.."grid/parallel_grid_layout_lv".. i .."_p"..ProcRank()..".ugx", offsetForGridHierarchy)


--------------------------
-- discretization setup --
--------------------------
-- scaled coefficient functions
function scaled_conductivity_out(z,r,t)
	local s = 2*math.pi*r * sigma_out
	return  s, 0,
			0, s
end
function scaled_conductivity_in(z,r,t)
	local s = 2*math.pi*r * sigma_in
	return  s, 0,
			0, s
end
function myelin_capacitance(z,r,t)
	return 2*math.pi*r * myelin_cap
end
function myelin_capacitance_neg(z,r,t)
	return - 2*math.pi*r * myelin_cap
end
function ranvier_capacitance(z,r,t)
	return 2*math.pi*r * ranvier_cap
end
function ranvier_capacitance_neg(z,r,t)
	return - 2*math.pi*r * ranvier_cap
end
function rotSym_scale(z,r,t)
	return 2*math.pi*r
end

-- outer potential
elemDiscPhiOut = ConvectionDiffusion("phiOut", "out", "fv1")
elemDiscPhiOut:set_diffusion("scaled_conductivity_out")
elemDiscPhiOut:set_mass_scale(0.0)


-- inner potential
elemDiscPhiIn = ConvectionDiffusion("phiIn", "in", "fv1")
elemDiscPhiIn:set_diffusion("scaled_conductivity_in")
elemDiscPhiIn:set_mass_scale(0.0)


-- gatings --
function zeroDeriv()
	return 0.0
end

-- n particle
elemDiscN = ConvectionDiffusion("n", ranvier, "fv1")

function rhsN(n, phiOut, phiIn)
	-- bit tricky: the result is computed in such a way
	-- that using it as (explicit) rhs in 'dn/dt = f' will make the solution
	-- of that ODE with the expl. Euler method
	-- the analytical solution for one time step (with constant Vm) 
	local vm = phiIn - phiOut
	local term = (n_inf(vm) - n)/dt
	return -term*(1.0-math.exp(-dt/tau_n(vm)))
end

rhsNData = LuaUserFunctionNumber("rhsN", 3)
rhsNData:set_input(0, elemDiscN:value())
rhsNData:set_input(1, elemDiscPhiOut:value())
rhsNData:set_input(2, elemDiscPhiIn:value())
rhsNData:set_deriv(0, "zeroDeriv")
rhsNData:set_deriv(1, "zeroDeriv")
rhsNData:set_deriv(2, "zeroDeriv")
elemDiscN:set_reaction(rhsNData)

-- m particle
elemDiscM = ConvectionDiffusion("m", ranvier, "fv1")

function rhsM(m, phiOut, phiIn)
	local vm = phiIn - phiOut
	local term = (m_inf(vm) - m)/dt
	return -term*(1.0 - math.exp(-dt/tau_m(vm))) 
end

rhsMData = LuaUserFunctionNumber("rhsM", 3)
rhsMData:set_input(0, elemDiscM:value())
rhsMData:set_input(1, elemDiscPhiOut:value())
rhsMData:set_input(2, elemDiscPhiIn:value())
rhsMData:set_deriv(0, "zeroDeriv")
rhsMData:set_deriv(1, "zeroDeriv")
rhsMData:set_deriv(2, "zeroDeriv")
elemDiscM:set_reaction(rhsMData)

-- h particle
elemDiscH = ConvectionDiffusion("h", ranvier, "fv1")

function rhsH(h, phiOut, phiIn)
	local vm = phiIn - phiOut
	local term = (h_inf(vm) - h)/dt
	return -term*(1.0-math.exp(-dt/tau_h(vm)))
end

rhsHData = LuaUserFunctionNumber("rhsH", 3)
rhsHData:set_input(0, elemDiscH:value())
rhsHData:set_input(1, elemDiscPhiOut:value())
rhsHData:set_input(2, elemDiscPhiIn:value())
rhsHData:set_deriv(0, "zeroDeriv")
rhsHData:set_deriv(1, "zeroDeriv")
rhsHData:set_deriv(2, "zeroDeriv")
elemDiscH:set_reaction(rhsHData)


-- Hodgkin-Huxley current
function hhCurrent(phiOut, phiIn, n, m, h, z, r, t, si)
	local potassiumFlux = g_K * n*n*n*n * (phiIn - phiOut - E_K)
	local sodiumFlux = g_Na * m*m*m*h * (phiIn - phiOut - E_Na)
	local leakageFlux = g_L * (phiIn - phiOut - E_L)
	return 2*math.pi*r * (potassiumFlux + sodiumFlux + leakageFlux)
end

hhCurrentFunction = LuaUserFunctionNumber("hhCurrent", 5, true)
hhCurrentFunction:set_input(0, elemDiscPhiOut:value())
hhCurrentFunction:set_input(1, elemDiscPhiIn:value())
hhCurrentFunction:set_input(2, elemDiscN:value())
hhCurrentFunction:set_input(3, elemDiscM:value())
hhCurrentFunction:set_input(4, elemDiscH:value())
hhCurrentFunction:set_deriv(0, "zeroDeriv")
hhCurrentFunction:set_deriv(1, "zeroDeriv")
hhCurrentFunction:set_deriv(2, "zeroDeriv")
hhCurrentFunction:set_deriv(3, "zeroDeriv")
hhCurrentFunction:set_deriv(4, "zeroDeriv")

linkerOut = ScaleAddLinkerNumber()
linkerOut:add(-1.0, hhCurrentFunction)
elemDiscHHOut = ConvectionDiffusionFV1("phiOut", ranvier)
elemDiscHHOut:set_mass_scale(0.0)
elemDiscHHOut:set_reaction(linkerOut)

linkerIn = ScaleAddLinkerNumber()
linkerIn:add(1.0, hhCurrentFunction)
elemDiscHHIn = ConvectionDiffusionFV1("phiIn", ranvier)
elemDiscHHIn:set_mass_scale(0.0)
elemDiscHHIn:set_reaction(linkerIn)


-- capacitive current
myelinCapData = LuaUserNumber("myelin_capacitance")
myelinCapDataNeg = LuaUserNumber("myelin_capacitance_neg")
ranvierCapData = LuaUserNumber("ranvier_capacitance")
ranvierCapDataNeg = LuaUserNumber("ranvier_capacitance_neg")

elemDiscCapMyelinOut = ConvectionDiffusion("phiOut", "myelin", "fv1")
elemDiscCapMyelinOut:set_mass_scale(myelinCapData)
linkerCapOutMyelin = ScaleAddLinkerNumber()
linkerCapOutMyelin:add(myelinCapDataNeg, elemDiscPhiIn:value())
elemDiscCapMyelinOut:set_mass(linkerCapOutMyelin)

elemDiscCapRanvierOut = ConvectionDiffusion("phiOut", ranvier, "fv1")
elemDiscCapRanvierOut:set_mass_scale(ranvierCapData)
linkerCapOutRanvier = ScaleAddLinkerNumber()
linkerCapOutRanvier:add(ranvierCapDataNeg, elemDiscPhiIn:value())
elemDiscCapRanvierOut:set_mass(linkerCapOutRanvier)

elemDiscCapMyelinIn = ConvectionDiffusion("phiIn", "myelin", "fv1")
elemDiscCapMyelinIn:set_mass_scale(myelinCapData)
linkerCapInMyelin = ScaleAddLinkerNumber()
linkerCapInMyelin:add(myelinCapDataNeg, elemDiscPhiOut:value())
elemDiscCapMyelinIn:set_mass(linkerCapInMyelin)

elemDiscCapRanvierIn = ConvectionDiffusion("phiIn", ranvier, "fv1")
elemDiscCapRanvierIn:set_mass_scale(ranvierCapData)
linkerCapInRanvier = ScaleAddLinkerNumber()
linkerCapInRanvier:add(ranvierCapDataNeg, elemDiscPhiOut:value())
elemDiscCapRanvierIn:set_mass(linkerCapInRanvier)


-- injection current
elemDiscInjection = UserFluxBoundaryFV1("phiIn, phiOut", injection)  -- function order: to, from
elemDiscInjection:set_flux_function("electrodeCurrent")


-- Dirichlet-0 boundary for outer potential
dirichletBndPhiOut = DirichletBoundary(true)
dirichletBndPhiOut:add(0.0, "phiOut", "bnd")


-- domain discretization
domainDiscImpl = DomainDiscretization(approxSpace)
domainDiscImpl:add(elemDiscPhiOut)
domainDiscImpl:add(elemDiscPhiIn)
domainDiscImpl:add(elemDiscCapMyelinOut)
domainDiscImpl:add(elemDiscCapRanvierOut)
domainDiscImpl:add(elemDiscCapMyelinIn)
domainDiscImpl:add(elemDiscCapRanvierIn)
domainDiscImpl:add(dirichletBndPhiOut)

domainDiscExpl = DomainDiscretization(approxSpace)
domainDiscExpl:add(elemDiscN)
domainDiscExpl:add(elemDiscM)
domainDiscExpl:add(elemDiscH)
domainDiscExpl:add(elemDiscHHOut)
domainDiscExpl:add(elemDiscHHIn)
domainDiscExpl:add(elemDiscInjection)  -- needs to be explicit, otherwise not assembled


-- time discretization
timeDiscImpl = ThetaTimeStep(domainDiscImpl)
timeDiscImpl:set_theta(1.0)

timeDiscExpl = ThetaTimeStep(domainDiscExpl)
timeDiscExpl:set_theta(0.0)

timeDisc = CompositeTimeDiscretization()
timeDisc:add_time_disc(timeDiscImpl)
timeDisc:add_time_disc(timeDiscExpl)


-------------------
-- algebra setup --
-------------------
-- linear operator
op = AssembledLinearOperator(timeDisc)

-- debug writer
dbgWriter = GridFunctionDebugWriter(approxSpace)
dbgWriter:set_vtk_output(false)
dbgWriter:set_base_dir(outDir)

convCheck = ConvCheck()
convCheck:set_maximum_steps(100)
convCheck:set_minimum_defect(1e-30)
convCheck:set_reduction(1e-08)
convCheck:set_verbose(true)
--[[
convCheck = CompositeConvCheck(approxSpace, 100, 1e-30, 1e-08)
convCheck:set_group_check("phiOut, phiIn", 1e-30, 1e-08)
convCheck:set_group_check("n, m, h", 1e-30, 1e-08)
convCheck:set_verbose(true)
convCheck:set_time_measurement(true)
--]]

-- solver setup --
solver = CG()--LinearSolver()

gmg = GeometricMultiGrid(approxSpace)
gmg:set_discretization(timeDisc)
gmg:set_base_level(0)
gmg:set_gathered_base_solver_if_ambiguous(true)
--gmg:set_debug(dbgWriter)

base = SuperLU()
gmg:set_base_solver(base)

smoother = ILU()--SymmetricGaussSeidel()
smoother:enable_consistent_interfaces(true)
smoother:set_sort(true) -- important since natural UG sorting is somewhat awkward
gmg:set_smoother(smoother)
gmg:set_smooth_on_surface_rim(true)
gmg:set_cycle_type(1)
gmg:set_num_presmooth(3)
gmg:set_num_postsmooth(3)
gmg:set_rap(true)

solver:set_preconditioner(gmg)
solver:set_convergence_check(convCheck)
--solver:set_debug(dbgWriter)


-------------------
-- Time stepping --
-------------------
-- initial values
time = 0
step = 0
u = GridFunction(approxSpace)
u:set(0.0)

InterpolateInner(phi_eq_out, u,"phiOut", volOut, time)
InterpolateInner(phi_eq_in, u, "phiIn", volIn, time)

n_eq = n_inf(phi_eq_in-phi_eq_out)
m_eq = m_inf(phi_eq_in-phi_eq_out)
h_eq = h_inf(phi_eq_in-phi_eq_out)
InterpolateInner(n_eq, u, "n", ranvier, time)
InterpolateInner(m_eq, u, "m", ranvier, time)
InterpolateInner(h_eq, u, "h", ranvier, time)


-- write start solutions
if generateVTKoutput then
	out = VTKOutput()
	out:print(outDir.."vtk/solution", u, step, time)
end

-- measure initial potential
take_measurement(u, time, ranvier, "phiIn, phiOut", outDir.."meas/meas")
take_measurement(u, time, "meas", "phiOut", outDir.."meas/bulk")


-- create new grid functions for old value and rhs
uOld = u:clone()
b = GridFunction(approxSpace)

-- store grid function in vector of old solutions
solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, time)


-- start the time stepping
for step = 1, numTimeSteps do
	print("++++++ TIMESTEP " .. step .. " BEGIN ++++++")

	-- prepare step
	timeDisc:prepare_step(solTimeSeries, dt)

	if step == 1 then
		-- assemble inner Vm problem (matrix and rhs)
		timeDisc:assemble_linear(op, b)
		
		-- init solver (invert matrix - the actual solving)
		solver:init(op)
	else
		-- assemble inner Vm problem (only rhs, matrix is const)
		timeDisc:assemble_rhs(b)
	end

	if not solver:apply(u, b) then
		print("Error: Solver for membrane potential problem did not converge.")
		if generateVTKoutput then
			out:write_time_pvd(outDir.."vtk/solution", u)
		end
		return
	end
			
	-- update to new time
	time = time + dt
		
	-- plot solution
	if generateVTKoutput and step % modulo == 0 then
		out:print(outDir.."vtk/solution", u, step/modulo, time)
	end
	
	-- measure potential
	take_measurement(u, time, ranvier, "phiIn, phiOut", outDir.."meas/meas")
	take_measurement(u, time, "meas", "phiOut", outDir.."meas/bulk")

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

