--------------------------------------------------------------------------------
-- 1D cable equation on axon with Ranvier nodes.                              --
--                                                                            --
-- Author: Markus Breit                                                       --
-- Date:   2018-09-04                                                         --
--------------------------------------------------------------------------------

ug_load_script("ug_util.lua")

AssertPluginsLoaded({"cable_neuron"})

-- choose algebra
InitUG(3, AlgebraType("CPU", 1))


------------------------------------
-- process command line arguments --
------------------------------------
-- grid
grid = util.GetParam("-grid", "grids/axon5_1d.ugx")
nRanvier = util.GetParam("-nRanv", 5)

-- time stepping
numTimeSteps =  util.GetParamNumber("-nSteps", 500, "number of timesteps")
dt = util.GetParamNumber("-dt", 2e-5, "time step in seconds")

-- spatial discretization
numRefs = util.GetParamNumber("-numRefs", 4, "number of refinements")

-- choose outfile directory
outDir = util.GetParam("-outName", "test/neuron")
outDir = outDir.."/"

-- generate vtk output? (every modulo-th step)
generateVTKoutput = util.HasParamOption("-vtk")
modulo = util.GetParamNumber("-modulo", 1) 


---------------------------------------
-- problem constants and stimulation --
---------------------------------------
-- resistivity (in units of (V s m) / C)
spec_res = 0.407224494  -- deduced from PNP model

-- membrane specific capacitance
ranvier_cap = 3.6120456108572e-3  -- 3.54167512e-03  -- in C / (V m^2)
myelin_cap = 2.4101772506797e-4  -- 10 layers with 2 membranes each

-- channel conductances
g_Na = 1.2e+03  -- in C / (V s m^2)
g_K = 3.6e+02  -- in C / (V s m^2)
g_L = 3.0e+00  -- in C / (V s m^2)

-- equilibrium potentials
phi_eq_in = -0.070  -- in V
E_Na = 0.0640  -- in V
E_K = -0.0940  -- in V
leak_K = 0.03092
leak_Na = -0.00293
E_L = phi_eq_in + (leak_K+leak_Na)/g_L

-- injection current density
injCurrentDensity = 3.0  -- in C / (s m^2)
injectionStart = 0  -- in s
injectionEnd = 1  -- in s


-------------------------------
-- approximation space setup --
-------------------------------
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

-- name some useful subset collections
ranvier = "ranvier1"
for i = 2,nRanvier do
	ranvier = ranvier .. ", ranvier" .. i
end
memAll = "myelin, " .. ranvier
injection = "ranvier1"

-- create approximation space
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("v", "Lagrange", 1)
approxSpace:init_top_surface()
approxSpace:init_levels()

OrderCuthillMcKee(approxSpace, true)

print()
print(dom:domain_info():to_string())

--SaveGridHierarchyTransformed(dom:grid(), dom:subset_handler(), outDir.."grid/refined_grid_hierarchy" .."_p" .. ProcRank() .. ".ugx", 2e-6)
--SaveParallelGridLayout(dom:grid(), fileName.."grid/parallel_grid_layout_lv".. i .."_p"..ProcRank()..".ugx", offsetForGridHierarchy)


--------------------------
-- discretization setup --
--------------------------
-- cable equation
cableMyelin = CableEquation("myelin", false)
cableMyelin:set_spec_cap(myelin_cap)
cableMyelin:set_spec_res(spec_res)
cableMyelin:set_rev_pot_k(E_K)
cableMyelin:set_rev_pot_na(E_Na)

cableRanvier = CableEquation(ranvier, false)
cableRanvier:set_spec_cap(ranvier_cap)
cableRanvier:set_spec_res(spec_res)
cableRanvier:set_rev_pot_k(E_K)
cableRanvier:set_rev_pot_na(E_Na)


-- Hodgkin and Huxley channels
hh = ChannelHH("v", ranvier)
hh:set_conductances(g_K, g_Na, ranvier)
hh:enable_temperature_dependency(false)
cableRanvier:add(hh)


-- leakage
leak = ChannelLeak("v", ranvier)
leak:set_cond(g_L, ranvier)
leak:set_rev_pot(E_L, ranvier)
cableRanvier:add(leak)


-- electrode stimulation (into subset 1, aka "ranvier1")
cableRanvier:set_influx_subset("ranvier1", injCurrentDensity, injectionEnd - injectionStart, injectionStart)


-- domain discretization
domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(cableMyelin)
domainDisc:add(cableRanvier)

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
dbgWriter:set_vtk_output(true)

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

Interpolate(phi_eq_in, u, "v")


-- write start solutions
if generateVTKoutput then
	out = VTKOutput()
	out:print(outDir.."vtk/solution", u, step, time)
end

-- measure initial potential
take_measurement(u, time, ranvier, "v", outDir.."meas/meas")


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
		exit()
	end
			
	-- update to new time
	time = time + dt
		
	-- plot solution
	if generateVTKoutput and step % modulo == 0 then
		out:print(outDir.."vtk/solution", u, step/modulo, time)
	end
	
	-- measure potential
	take_measurement(u, time, ranvier, "v", outDir.."meas/meas")

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

