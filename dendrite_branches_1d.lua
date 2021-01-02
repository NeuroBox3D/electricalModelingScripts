--------------------------------------------------------------------------------
-- 3D cable equation on branching dendrite.                                   --
-- The outer membrane potential is supposed to be a constant zero.            --
--                                                                            --
-- Author: Markus Breit                                                       --
-- Date:   2018-12-28                                                         --
--------------------------------------------------------------------------------

ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")

AssertPluginsLoaded({"cable_neuron"})

-- choose algebra
InitUG(3, AlgebraType("CPU", 1))


------------------------------------
-- process command line arguments --
------------------------------------
-- grid
grid = util.GetParam("-grid", "grids/branches_1d.ugx")

-- time stepping
numTimeSteps =  util.GetParamNumber("-nSteps", 500, "number of timesteps")
dt = util.GetParamNumber("-dt", 2e-5, "time step in seconds")

-- spatial discretization
numRefs = util.GetParamNumber("-numRefs", 0, "number of refinements")

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
spec_res = 0.407224494  -- deduced from PNP model, in (V s m) / C

-- membrane specific capacitance
cap = 3.6120456108572e-03  -- match params from neuron_mb.lua, otherwise 1e-02  -- in C / (V m^2)

-- equilibrium potential
vm_eq = -0.07  -- in V

-- leakage
g_L = 3.0  -- in C / (V s m^2)

-- injection current density
injCurrentDensity = 5.0  -- in C / (s m^2)
injectionStart = 0  -- in s
injectionEnd = 0.005  -- in s


-------------------------------
-- approximation space setup --
-------------------------------
-- create domain
dom = Domain()
dom:create_additional_subset_handler("projSH")
LoadDomain(dom, grid)

-- check subsets
neededSubsets = {"in", "syn", "meas"}
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

-- create approximation space
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("v", "Lagrange", 1)
approxSpace:init_top_surface()
approxSpace:init_levels()

OrderCuthillMcKee(approxSpace, true)

print()
print(dom:domain_info():to_string())

--SaveGridHierarchyTransformed(dom:grid(), dom:subset_handler(), outDir.."grid/refined_grid_hierarchy" .."_p" .. ProcRank() .. ".ugx", 2.0)
--SaveParallelGridLayout(dom:grid(), fileName.."grid/parallel_grid_layout_lv".. i .."_p"..ProcRank()..".ugx", offsetForGridHierarchy)


--------------------------
-- discretization setup --
--------------------------
-- cable equation
cableDisc = CableEquation("in, syn, meas", false)
cableDisc:set_spec_cap(cap)
cableDisc:set_spec_res(spec_res)

-- leakage
leak = ChannelLeak("v", "in, syn, meas")
leak:set_cond(g_L, "in, syn, meas")
leak:set_rev_pot(vm_eq, "in, syn, meas")
cableDisc:add(leak)

-- electrode stimulation (into subset 1, aka "syn")
cableDisc:set_influx_subset("syn", injCurrentDensity, injectionEnd - injectionStart, injectionStart)


-- domain discretization
domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(cableDisc)

-- time discretization
timeDisc = ThetaTimeStep(domainDisc)
timeDisc:set_theta(1.0)


-- create linear operator
op = AssembledLinearOperator(timeDisc)


-------------------
-- algebra setup --
-------------------
-- debug writer
dbgWriter = GridFunctionDebugWriter(approxSpace)
dbgWriter:set_vtk_output(false)
dbgWriter:set_base_dir(outDir)

-- linear solver --
convCheck = CompositeConvCheck(approxSpace, 20, 2e-26, 1e-08)
convCheck:set_component_check("v", 1e-21, 1e-12)
convCheck:set_verbose(true)

ilu = ILU()
solver = LinearSolver()
solver:set_preconditioner(ilu)
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

Interpolate(vm_eq, u, "v")


-- write start solutions
if generateVTKoutput then
	out = VTKOutput()
	out:print(outDir.."vtk/solution", u, step, time)
end

-- measure initial potential
take_measurement(u, time, "meas", "v", outDir.."meas/meas")


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
	take_measurement(u, time, "meas", "v", outDir.."meas/meas")

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

