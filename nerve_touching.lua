--------------------------------------------------------------------------------
-- 3D cable equation on bundle of 7 touching nerve fibers with Ranvier nodes. --
-- The fibres are touching one another.                                       --
-- Parallel execution is required (~100-1000 procs).                          --
--                                                                            --
-- Author: Markus Breit                                                       --
-- Date:   2018-09-14                                                         --
--------------------------------------------------------------------------------

ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")

AssertPluginsLoaded({"ConvectionDiffusion", "neuro_collection",
					 "neuro_collection", "Parmetis"}) 

-- choose algebra
InitUG(3, AlgebraType("CPU", 1))


------------------------------------
-- process command line arguments --
------------------------------------
-- grid
grid = util.GetParam("-grid", "grids/nerve00.ugx")
nRanvier = util.GetParam("-nRanv", 5)

-- time stepping
numTimeSteps =  util.GetParamNumber("-nSteps", 500, "number of timesteps")
dt = util.GetParamNumber("-dt", 0.02, "time step in milliseconds")

-- spatial discretization
numRefs = util.GetParamNumber("-numRefs", 2, "number of refinements")
numAnisoRefs = util.GetParamNumber("-numAnisoRefs", 5, "number of anisotropic refinements")

-- stimulation of central or surrounding (default) neurons
centralStim = util.HasParamOption("-centralStim")

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
sigma_out = 2.4556479647688e3  -- deduced from PNP model, in fC / (ms mV um)

-- membrane specific capacitance
ranvier_cap = 3.6120456108572e-3  -- match params from neuron_mb.lua, in fC / (mV um^2)
myelin_cap = 2.4101772506797e-4  -- 10 layers with 2 membranes each

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
E_L = (phi_eq_in - phi_eq_out) + (leak_K+leak_Na)/g_L

-- injection current density
injCurrentDensity = 3.0  -- in fC / (ms um^2)


-------------------------------
-- approximation space setup --
-------------------------------
-- create domain
dom = Domain()
dom:create_additional_subset_handler("projSH")
LoadDomain(dom, grid)

-- check subsets
neededSubsets = {"bnd"}
for j = 0, 6 do neededSubsets[1+j] = "in" .. j end
neededSubsets[8] = "out"
for j = 0, 6 do neededSubsets[9+j] = "myelin" .. j end
for i = 1, nRanvier do
	for j = 0, 6 do
		neededSubsets[16+7*(i-1)+j] = "ranvier" .. i .. j
	end
end
neededSubsets[16+7*nRanvier] = "bnd"
ug_assert(util.CheckSubsets(dom, neededSubsets), "Something wrong with required subsets.")

-- distribute
balancer.partitioner = "parmetis"
balancer.staticProcHierarchy = true
balancer.firstDistLvl = (numAnisoRefs+numRefs+1-math.abs(numAnisoRefs+numRefs-1))/2
balancer.redistSteps = 2
balancer.firstDistProcs = 120
balancer.redistProcs = 16
balancer.parallelElementThreshold = 24
balancer.ParseParameters()

loadBalancer = balancer.CreateLoadBalancer(dom)
if loadBalancer ~= nil then
	if balancer.partitioner == "parmetis" then
		-- new Parmetis implementation
		mu = ManifoldUnificator(dom)
		for j = 0, 6 do
			mu:add_protectable_subsets("myelin"..j)
			for i = 1, nRanvier do
				mu:add_protectable_subsets("ranvier" .. i .. j)
			end
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
	loadBalancer:print_quality_records()
end

-- name some useful subset collections
inner = {}
myelin = {}
myelinExtra = {}
ranvier = {}
ranvierExtra = {}
memAll = {}
volIn = {}
volOut = "out, bnd"
for j = 0, 6 do
	inner[j+1] = "in" .. j
	myelin[j+1] = "myelin" .. j
	ranvier[j+1] = "ranvier" .. 1 .. j	
	for i = 2, nRanvier do
		ranvier[j+1] = ranvier[j+1] .. ", ranvier" .. i .. j
	end
	
	myelinExtra[j+1] = ""
	ranvierExtra[j+1] = ""
	if j == 0 then
		for k = 1, 6 do
			myelinExtra[j+1] = myelinExtra[j+1] .. ", myelin" .. j .. k
			for i = 1, nRanvier do
				ranvierExtra[j+1] = ranvierExtra[j+1] .. ", ranvier" .. i .. j .. k
			end
		end
	else
		myelinExtra[j+1] = myelinExtra[j+1] .. ", myelin" .. 0 .. j
		for i = 1, nRanvier do
			ranvierExtra[j+1] = ranvierExtra[j+1] .. ", ranvier" .. i .. 0 .. j
		end
		local last = (j+4) % 6 + 1
		local next = j%6 + 1
		myelinExtra[j+1] = myelinExtra[j+1] .. ", myelin" .. last .. j
		myelinExtra[j+1] = myelinExtra[j+1] .. ", myelin" .. j .. next
		for i = 1, nRanvier do
			ranvierExtra[j+1] = ranvierExtra[j+1] .. ", ranvier" .. i .. last .. j
			ranvierExtra[j+1] = ranvierExtra[j+1] .. ", ranvier" .. i .. j .. next
		end
	end
	
	memAll[j+1] = myelin[j+1] .. myelinExtra[j+1] .. ", " .. ranvier[j+1] .. ranvierExtra[j+1]
	volIn[j+1] = inner[j+1] .. ", " .. memAll[j+1]
	volOut = volOut .. ", " .. memAll[j+1]
end

injection = "ranvier10"

--[[
for j = 0, 6 do
	print("myelin[" .. j .. "] = " .. myelin[j+1])
	print("ranvier[" .. j .. "] = " .. ranvier[j+1])
	print("myelinExtra[" .. j .. "] = " .. myelinExtra[j+1])
	print("ranvierExtra[" .. j .. "] = " .. ranvierExtra[j+1])
	print("memAll[" .. j .. "] = " .. memAll[j+1])
	print("volIn[" .. j .. "] = " .. volIn[j+1])
end
print("volOut = " .. volOut)
--]]

meas = ""
for i = 1, nRanvier do
	meas = meas .. ", ranvier" .. i .. 0
end
meas = string.sub(meas, 3)


-- create approximation space
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("phiOut", "Lagrange", 1, volOut)
for j = 0, 6 do
	approxSpace:add_fct("phiIn"..j, "Lagrange", 1, volIn[j+1])
	approxSpace:add_fct("n"..j, "Lagrange", 1, ranvier[j+1] .. ranvierExtra[j+1])
	approxSpace:add_fct("m"..j, "Lagrange", 1, ranvier[j+1] .. ranvierExtra[j+1])
	approxSpace:add_fct("h"..j, "Lagrange", 1, ranvier[j+1] .. ranvierExtra[j+1])
end
approxSpace:init_top_surface()
approxSpace:init_levels()

-- anisotropic refinements
refiner = HangingNodeDomainRefiner(dom)
for j = 0, 6 do
	local rma = ExtensionRefMarkAdjuster(dom, {1,0,0}, "in" .. j)
	add_extension_ref_mark_adjuster(refiner, rma)
end
rmaOut = ExtensionRefMarkAdjuster(dom, {1,0,0}, "out")
add_extension_ref_mark_adjuster(refiner, rmaOut)
for i = 1, numAnisoRefs do
	MarkGlobal(refiner, dom)
	refiner:refine()
	balancer.Rebalance(dom, loadBalancer)
end

-- isotropic refinements
refiner = GlobalDomainRefiner(dom)
for i = 1, numRefs do
	refiner:refine()
	balancer.Rebalance(dom, loadBalancer)
end

print()
print(dom:domain_info():to_string())

--SaveGridHierarchyTransformed(dom:grid(), dom:subset_handler(), outDir.."grid/refined_grid_hierarchy" .."_p" .. ProcRank() .. ".ugx", 3.0)
--SaveParallelGridLayout(dom:grid(), fileName.."grid/parallel_grid_layout_lv".. i .."_p"..ProcRank()..".ugx", 3.0)


--------------------------
-- discretization setup --
--------------------------
-- outer potential
elemDiscPhiOut = ConvectionDiffusion("phiOut", "out", "fv1")
elemDiscPhiOut:set_diffusion(sigma_out)
elemDiscPhiOut:set_mass_scale(0.0)

-- inner potential
elemDiscPhiIn = {}
hh = {}
elemDiscHH = {}
leakage = {}
elemDiscLeakage = {}
elemDiscCapMyelinOut = {}
linkerCapOutMyelin = {}
elemDiscCapRanvierOut = {}
linkerCapOutRanvier = {}
elemDiscCapMyelinIn = {}
linkerCapInMyelin = {}
elemDiscCapRanvierIn = {}
linkerCapInRanvier = {}

for j = 0, 6 do
	elemDiscPhiIn[j+1] = ConvectionDiffusion("phiIn"..j, inner[j+1], "fv1")
	elemDiscPhiIn[j+1]:set_diffusion(sigma_in)
	elemDiscPhiIn[j+1]:set_mass_scale(0.0)

	-- Hodgkin-Huxley channels
	-- let's pretend we had a constant density of HH K and Na channels of 1/um^2
	-- (since we do not know how membrane conductance is distributed between density and single-channel conductance)
	hh[j+1] = HH("phiIn" .. j .. ", phiOut, n" .. j .. ", m" .. j .. ", h" .. j, ranvier[j+1])
	hh[j+1]:set_conductances(1e3*g_K, 1e3*g_Na)
	hh[j+1]:set_reversal_potentials(1e-3*E_K, 1e-3*E_Na)
	hh[j+1]:set_reference_time(1e-3)
	hh[j+1]:set_scale_inputs({1e-3, 1e-3, 1.0, 1.0, 1.0})
	hh[j+1]:use_exact_gating_mode(dt)
	elemDiscHH[j+1] = MembraneTransportFV1(ranvier[j+1], hh[j+1])
	elemDiscHH[j+1]:set_density_function(1.0)
	
	-- leakage
	leakage[j+1] = OhmicLeakage("phiIn" .. j .. ", phiOut")
	leakage[j+1]:set_conductance(1e3*g_L)
	leakage[j+1]:set_reversal_potential(1e-3*E_L)
	leakage[j+1]:set_scale_inputs({1e-3, 1e-3})
	elemDiscLeakage[j+1] = MembraneTransportFV1(ranvier[j+1], leakage[j+1])
	elemDiscLeakage[j+1]:set_density_function(1.0)
	
	-- capacitive current
	elemDiscCapMyelinOut[j+1] = ConvectionDiffusionFV1("phiOut", myelin[j+1])
	elemDiscCapMyelinOut[j+1]:set_mass_scale(myelin_cap)
	linkerCapOutMyelin[j+1] = ScaleAddLinkerNumber()
	linkerCapOutMyelin[j+1]:add(-myelin_cap, elemDiscPhiIn[j+1]:value())
	elemDiscCapMyelinOut[j+1]:set_mass(linkerCapOutMyelin[j+1])
	
	elemDiscCapRanvierOut[j+1] = ConvectionDiffusionFV1("phiOut", ranvier[j+1])
	elemDiscCapRanvierOut[j+1]:set_mass_scale(ranvier_cap)
	linkerCapOutRanvier[j+1] = ScaleAddLinkerNumber()
	linkerCapOutRanvier[j+1]:add(-ranvier_cap, elemDiscPhiIn[j+1]:value())
	elemDiscCapRanvierOut[j+1]:set_mass(linkerCapOutRanvier[j+1])
	
	elemDiscCapMyelinIn[j+1] = ConvectionDiffusionFV1("phiIn"..j, myelin[j+1])
	elemDiscCapMyelinIn[j+1]:set_mass_scale(myelin_cap)
	linkerCapInMyelin[j+1] = ScaleAddLinkerNumber()
	linkerCapInMyelin[j+1]:add(-myelin_cap, elemDiscPhiOut:value())
	elemDiscCapMyelinIn[j+1]:set_mass(linkerCapInMyelin[j+1])
	
	elemDiscCapRanvierIn[j+1] = ConvectionDiffusionFV1("phiIn"..j, ranvier[j+1])
	elemDiscCapRanvierIn[j+1]:set_mass_scale(ranvier_cap)
	linkerCapInRanvier[j+1] = ScaleAddLinkerNumber()
	linkerCapInRanvier[j+1]:add(-ranvier_cap, elemDiscPhiOut:value())
	elemDiscCapRanvierIn[j+1]:set_mass(linkerCapInRanvier[j+1])

end

-- injection current
if centralStim then
	elemDiscInjection = UserFluxBoundaryFV1("phiIn0, phiOut", injection)  -- function order: to, from
	elemDiscInjection:set_flux_function(injCurrentDensity)
else
	elemDiscInjection = {}
	for j = 1, 6 do
		elemDiscInjection[j] = UserFluxBoundaryFV1("phiIn"..j..", phiOut", "ranvier1"..j)  -- function order: to, from
		elemDiscInjection[j]:set_flux_function(injCurrentDensity)
	end
end

-- Dirichlet-0 boundary for outer potential
dirichletBndPhiOut = DirichletBoundary(true)
dirichletBndPhiOut:add(0.0, "phiOut", "bnd")


-- domain discretization
domainDiscImpl = DomainDiscretization(approxSpace)
domainDiscImpl:add(elemDiscPhiOut)
for j = 0, 6 do
	domainDiscImpl:add(elemDiscPhiIn[j+1])
	domainDiscImpl:add(elemDiscCapMyelinOut[j+1])
	domainDiscImpl:add(elemDiscCapRanvierOut[j+1])
	domainDiscImpl:add(elemDiscCapMyelinIn[j+1])
	domainDiscImpl:add(elemDiscCapRanvierIn[j+1])
end
domainDiscImpl:add(dirichletBndPhiOut)

domainDiscExpl = DomainDiscretization(approxSpace)
for j = 0, 6 do
	domainDiscExpl:add(hh[j+1])
	domainDiscExpl:add(elemDiscHH[j+1])
	domainDiscExpl:add(elemDiscLeakage[j+1])
end

if centralStim then
	domainDiscExpl:add(elemDiscInjection)  -- needs to be explicit, otherwise not assembled
else
	for j = 1, 6 do
		domainDiscExpl:add(elemDiscInjection[j])
	end
end

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
gmg:set_base_level((numAnisoRefs+numRefs+1-math.abs(numAnisoRefs+numRefs-1))/2)
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

InterpolateInner(phi_eq_out, u, "phiOut", volOut, time)
for j = 0, 6 do
	InterpolateInner(phi_eq_in, u, "phiIn"..j, volIn[j+1], time)
end

n_eq = n_inf(phi_eq_in-phi_eq_out)
m_eq = m_inf(phi_eq_in-phi_eq_out)
h_eq = h_inf(phi_eq_in-phi_eq_out)
for j = 0, 6 do
	InterpolateInner(n_eq, u, "n"..j, ranvier[j+1]..ranvierExtra[j+1], time)
	InterpolateInner(m_eq, u, "m"..j, ranvier[j+1]..ranvierExtra[j+1], time)
	InterpolateInner(h_eq, u, "h"..j, ranvier[j+1]..ranvierExtra[j+1], time)
end

-- write start solutions
if generateVTKoutput then
	out = VTKOutput()
	out:print(outDir.."vtk/solution", u, step, time)
end

-- measure initial potential
take_measurement(u, time, meas, "phiIn0, phiOut", outDir.."meas/meas")


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
	take_measurement(u, time, meas, "phiIn0, phiOut", outDir.."meas/meas")

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

