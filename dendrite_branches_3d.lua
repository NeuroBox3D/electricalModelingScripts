--------------------------------------------------------------------------------
-- 3D cable equation on branching dendrite.                                   --
-- The outer membrane potential is supposed to be a constant zero.            --
--                                                                            --
-- Author: Markus Breit                                                       --
-- Date:   2018-12-28                                                         --
--------------------------------------------------------------------------------

ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")

AssertPluginsLoaded({"ConvectionDiffusion", "Parmetis"}) 

-- choose algebra
InitUG(3, AlgebraType("CPU", 1))


------------------------------------
-- process command line arguments --
------------------------------------
-- grid
grid = util.GetParam("-grid", "grids/branches_3d.ugx")

-- time stepping
numTimeSteps =  util.GetParamNumber("-nSteps", 500, "number of timesteps")
dt = util.GetParamNumber("-dt", 0.02, "time step in milliseconds")

-- spatial discretization
numRefs = util.GetParamNumber("-numRefs", 3, "number of refinements")
numAnisoRefs = util.GetParamNumber("-numAnisoRefs", 0, "number of refinements")

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
sigma = 2.4556479647688e3  -- deduced from PNP model, in fC / (ms mV um)

-- membrane specific capacitance
cap = 3.54167512e-03  -- match params from neuron_mb.lua, otherwise 1e-02  -- in fC / (mV um^2)

-- equilibrium potential
vm_eq = -70.0  -- in mV

-- leakage
g_L = 0.003  -- in fC / (mV ms um^2)

-- injection current density
injCurrentDensity = 5.0  -- in fC / (ms um^2)
injectionStart = 0  -- in ms
injectionEnd = 5  -- in ms

function electrodeCurrent(x, y, z, t, si)
	if t >= injectionStart and t <= injectionEnd then
		return injCurrentDensity
	end
	return 0.0
end


-------------------------------
-- approximation space setup --
-------------------------------
-- create domain
dom = Domain()
dom:create_additional_subset_handler("projSH")
LoadDomain(dom, grid)

-- check subsets
neededSubsets = {"in", "mem", "syn", "meas"}
ug_assert(util.CheckSubsets(dom, neededSubsets), "Something wrong with required subsets.")

-- distribute
balancer.partitioner = "parmetis"
balancer.staticProcHierarchy = true
balancer.firstDistLvl = -1
balancer.redistSteps = 0
balancer.parallelElementThreshold = 24
balancer.ParseParameters()

loadBalancer = balancer.CreateLoadBalancer(dom)
if loadBalancer ~= nil then
	balancer.Rebalance(dom, loadBalancer)
	loadBalancer:estimate_distribution_quality()
	loadBalancer:print_quality_records()
end


-- create approximation space
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("vm", "Lagrange", 1)
approxSpace:init_top_surface()
approxSpace:init_levels()

-- anisotropic refinements
--refiner = HangingNodeDomainRefiner(dom)
-- TODO

-- isotropic refinements
refiner = GlobalDomainRefiner(dom)
for i = 1, numRefs do
	refiner:refine()
end


print()
print(dom:domain_info():to_string())

SaveGridHierarchyTransformed(dom:grid(), dom:subset_handler(), outDir.."grid/refined_grid_hierarchy" .."_p" .. ProcRank() .. ".ugx", 2.0)
--SaveParallelGridLayout(dom:grid(), fileName.."grid/parallel_grid_layout_lv".. i .."_p"..ProcRank()..".ugx", offsetForGridHierarchy)


--------------------------
-- discretization setup --
--------------------------
-- potential
elemDiscVm = ConvectionDiffusionFV1("vm", "in")
elemDiscVm:set_diffusion(sigma)
elemDiscVm:set_mass_scale(0.0)

-- leakage
elemDiscLeakage = ConvectionDiffusionFV1("vm", "mem, syn, meas")
elemDiscLeakage:set_mass_scale(0.0)
elemDiscLeakage:set_reaction_rate(g_L)
elemDiscLeakage:set_source(g_L*vm_eq)

-- capacitive current
elemDiscCapacitive = ConvectionDiffusionFV1("vm", "mem, syn, meas")
elemDiscCapacitive:set_mass_scale(cap)

-- injection current
elemDiscInjection = UserFluxBoundaryFV1("vm", "syn")  -- function order: to, from
elemDiscInjection:set_flux_function("electrodeCurrent")


-- domain discretization
domainDiscImpl = DomainDiscretization(approxSpace)
domainDiscImpl:add(elemDiscVm)
domainDiscImpl:add(elemDiscCapacitive)
domainDiscImpl:add(elemDiscLeakage)

domainDiscExpl = DomainDiscretization(approxSpace)
domainDiscExpl:add(elemDiscInjection)  -- needs to be explicit, otherwise not assembled


-- time discretization
timeDiscImpl = ThetaTimeStep(domainDiscImpl)
timeDiscImpl:set_theta(1.0)

timeDiscExpl = ThetaTimeStep(domainDiscExpl)
timeDiscExpl:set_theta(0.0)

timeDisc = CompositeTimeDiscretization()
timeDisc:add_time_disc(timeDiscImpl)
timeDisc:add_time_disc(timeDiscExpl)


-- create linear operator
opLinear = AssembledLinearOperator(timeDisc)


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
convCheck:set_group_check("vm", 1e-30, 1e-08)
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
InterpolateInner(vm_eq, u, "vm", "in, mem, syn, meas", time)


-- write start solutions
if generateVTKoutput then
	out = VTKOutput()
	out:print(outDir.."vtk/solution", u, step, time)
end

-- measure initial potential
take_measurement(u, time, "meas", "vm", outDir.."meas/meas")


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
		out:write_time_pvd(outDir.."vtk/solution", u)
		exit()
	end
			
	-- update to new time
	time = time + dt
		
	-- plot solution
	if generateVTKoutput and step % modulo == 0 then
		out:print(outDir.."vtk/solution", u, step/modulo, time)
	end
	
	-- measure potential
	take_measurement(u, time, "meas", "vm", outDir.."meas/meas")

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

