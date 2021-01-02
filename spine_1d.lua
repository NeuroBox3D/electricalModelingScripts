--------------------------------------------------------------------------------
-- This script solves the cable equation on a 1d version of a reconstructed   --
-- spine (with extensions).                                                   --
--                                                                            --
-- Author: Markus Breit                                                       --
-- Date:   2019-01-17                                                         --
--------------------------------------------------------------------------------
ug_load_script("ug_util.lua")

-- init UG
InitUG(3, AlgebraType("CPU", 1))
AssertPluginsLoaded({"cable_neuron"})


--------------
-- settings --
--------------
-- grid
gridName = util.GetParam("-grid", "grids/spine1d.ugx")

-- refinement
numRefs = util.GetParamNumber("-numRefs", 0)

-- time stepping
dt = util.GetParamNumber("-dt", 2e-5)  -- in s
endTime = util.GetParamNumber("-endTime", 0.02)  -- in s
nSteps = util.GetParamNumber("-nSteps", endTime/dt)

-- with simulation of single ion concentrations?
withIons = util.HasParamOption("-ions")

-- specify "-verbose" to output linear solver convergence
verbose	= util.HasParamOption("-verbose")

-- vtk output?
generateVTKoutput = util.HasParamOption("-vtk")
pstep = util.GetParamNumber("-pstep", dt)

-- file handling
outDir = util.GetParam("-outName", "viet_1d")
outDir = outDir .. "/"


--------------------------
-- biological settings	--
--------------------------
-- membrane conductances (in units of S/m^2)
g_l_k = 0.5
g_l_na = 0.5

-- resistivity (in units of (V s m) / C)
spec_res = 0.407224494  -- deduced from PNP model

-- membrane specific capacitance (in units of C / (V m^2))
spec_cap = 3.585763867e-3  -- deduced from PNP model

-- reversal potentials (in units of V)
e_k  = -0.0940
e_na = 0.0640
e_ca = 0.14

-- equilibrium concentrations (in units of mM)
k_out  = 4.0
na_out = 145.0
ca_out = 2.0

k_in = 155.0
na_in = 12.0
ca_in = 5e-5

-- equilibrium potential (in units of V)
v_eq = -0.070

-- diffusion coefficients (in units of m^2/s)
diff_k 	= 1.96e-9
diff_na	= 1.33e-9
diff_ca = 2.2e-10

-- temperature in units of K
temp = 298.15


--------------------------------
-- create approximation space --
--------------------------------
dom = util.CreateDomain(gridName, numRefs, {"dend"})

approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("v", "Lagrange", 1)
if withIons then
	approxSpace:add_fct("k", "Lagrange", 1)
	approxSpace:add_fct("na", "Lagrange", 1)
	approxSpace:add_fct("ca", "Lagrange", 1)
end
approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_statistic()
OrderCuthillMcKee(approxSpace, true)


-- cable equation
CE = CableEquation("dend, spine, syn", withIons)
CE:set_spec_cap(spec_cap)
CE:set_spec_res(spec_res)
CE:set_rev_pot_k(e_k)
CE:set_rev_pot_na(e_na)
CE:set_rev_pot_ca(e_ca)
CE:set_k_out(k_out)
CE:set_na_out(na_out)
CE:set_ca_out(ca_out)
CE:set_diff_coeffs({diff_k, diff_na, diff_ca})
CE:set_temperature(temp)

-- leakage
if withIons then
	leakK = IonLeakage("k", "dend")
	leakK:set_ohmic(true)
	leakK:set_valency(1)
	leakK:set_cond(g_l_k)
	leakK:set_rev_pot(v_eq)
	CE:add(leakK)
	
	leakNa = IonLeakage("na", "dend")
	leakNa:set_ohmic(true)
	leakNa:set_valency(1)
	leakNa:set_cond(g_l_na)
	leakNa:set_rev_pot(v_eq)
	CE:add(leakNa)
else
	leakK = ChannelLeak("k", "dend")
	leakK:set_cond(g_l_k, "dend")
	leakK:set_rev_pot(v_eq, "dend")
	CE:add(leakK)
	
	leakNa = ChannelLeak("na", "dend")
	leakNa:set_cond(g_l_na, "dend")
	leakNa:set_rev_pot(v_eq, "dend")
	CE:add(leakNa)
end

-- synapse
local t_on = 0.0
local tau = 2.5e-3  -- in s
local period = 0.01
local tau_nmdar = 0.15  -- in s
function channelOpen(t)
	-- strong NMDAR
	--return math.exp(-t/tau_nmdar) -- between 1/10 and 1/100 max
	
	-- tetanic stimulation
	local t_p = t % period
	--return math.exp(-t_p / tau)
	
	-- AMPAR
	return (t - t_on) / tau * math.exp(-(t - t_on - tau) / tau)
end

FRT = 96485.0 / (8.31451 * 298.15)  -- in 1 / V
synArea = 1.27971e-07
areaFactor = 9.275e-14 / (math.pi * 0.5*(2.95945e-07 + 1.78e-07) * synArea)
vmAtSyn = v_eq
function currentDensityFunction(x, y, z, t, si)
	local exp = math.exp(FRT*vmAtSyn)
	local currentK = 96485.0 * diff_k / 1e-8 * FRT * vmAtSyn
		* (k_out - k_in * exp) / (1.0 - exp)
	local currentNa = 96485.0 * diff_na / 1e-8 * FRT * vmAtSyn
		* (na_out - na_in * exp) / (1.0 - exp)
	return - channelOpen(t) * 1e-4 * (currentK + currentNa) * areaFactor
end

if withIons then
	synK = IonLeakage("k", "syn")
	synK:set_valency(1)
	synK:set_perm(0.0)
	CE:add(synK)
	
	synNa = IonLeakage("na", "syn")
	synNa:set_valency(1)
	synNa:set_perm(0.0)
	CE:add(synNa)
else
	CE:set_influx_function("currentDensityFunction", "syn")
end

-- create domain discretization
domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(CE)

assTuner = domainDisc:ass_tuner()


-- create time discretization
timeDisc = ThetaTimeStep(domainDisc)
timeDisc:set_theta(1.0)


-- create operator from discretization
linOp = AssembledLinearOperator(timeDisc)

------------------
-- solver setup	--
------------------
-- debug writer
dbgWriter = GridFunctionDebugWriter(approxSpace)
dbgWriter:set_vtk_output(true)

-- linear solver --
convCheck = CompositeConvCheck(approxSpace, 1, 2e-26, 1e-08)
convCheck:set_component_check("v", 1e-21, 1e-12)
convCheck:set_verbose(verbose)

ilu = ILU()
solver = LinearSolver()
solver:set_preconditioner(ilu)
solver:set_convergence_check(convCheck)
--solver:set_debug(dbgWriter)


----------------------
-- time stepping	--
----------------------
time = 0.0

-- init solution
u = GridFunction(approxSpace)
b = GridFunction(approxSpace)
u:set(0.0)
Interpolate(v_eq, u, "v")
if withIons then
	Interpolate(k_in, u, "k")
	Interpolate(na_in, u, "na")
	Interpolate(ca_in, u, "ca")
end

synArea = Integral(1.0, u, "syn")
areaFactor = 9.275e-14 / (math.pi * 0.5*(2.95945e-07 + 1.78e-07) * synArea)
currentOutFile = assert(io.open(outDir.."meas/current.dat", "a"))

-- write start solution
if generateVTKoutput then 
	out = VTKOutput()
	out:print(outDir.."vtk/solution", u, 0, time)
end
measFcts = "v"
if withIons then
	measFcts = "v, k, na"
end
take_measurement(u, time, "syn", measFcts, outDir.."meas/sol")

-- store grid function in vector of  old solutions
uOld = u:clone()
solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, time)

curr_dt = dt
dtred = 2

lv = 0
maxLv = 10
cb_counter = {}
cb_counter[lv] = 0
while endTime-time > 0.001*curr_dt do
		-- setup time Disc for old solutions and timestep
	timeDisc:prepare_step(solTimeSeries, curr_dt)
	
	vmAtSyn = Integral(u, "v", "syn") / synArea
	currentOutFile:write(time, "\t", -currentDensityFunction(0,0,0,time,0) * math.pi * 0.5*(2.95945e-07 + 1.78e-07) * synArea, "\n")
	if withIons then
		synK:set_perm(channelOpen(time) * 1e-4 * diff_k / 1e-8 * areaFactor)
		synNa:set_perm(channelOpen(time) * 1e-4 * diff_na / 1e-8 * areaFactor)
	end
	
	-- reduce time step if cfl < curr_dt
	-- (this needs to be done AFTER prepare_step as channels are updated there)
	dtChanged = false
	cfl = CE:estimate_cfl_cond(solTimeSeries:latest())
	print("estimated CFL condition: dt < " .. cfl)
	while (curr_dt > cfl) do
		curr_dt = curr_dt/dtred
		
		if lv+1 > maxLv then
			print("Time step too small.")
			exit()
		end
		
		lv = lv + 1
		cb_counter[lv] = 0
		print("estimated CFL condition: dt < " .. cfl .. " - reducing time step to " .. curr_dt)
		dtChanged = true
	end
	
	-- increase time step if cfl > curr_dt / dtred (and if time is aligned with new bigger step size)
	while curr_dt*dtred < cfl and lv > 0 and cb_counter[lv] % (dtred) == 0 do
		curr_dt = curr_dt*dtred
		lv = lv - 1
		cb_counter[lv] = cb_counter[lv] + cb_counter[lv+1]/dtred
		cb_counter[lv+1] = 0
		print ("estimated CFL condition: dt < " .. cfl .. " - increasing time step to " .. curr_dt)
		dtChanged = true
	end
	
	print("++++++ POINT IN TIME " .. math.floor((time+curr_dt)/curr_dt+0.5)*curr_dt .. " BEGIN ++++++")
	
	-- prepare again with new time step size
	if dtChanged == true then 
		timeDisc:prepare_step(solTimeSeries, curr_dt)
	end

	-- assemble linear problem
	matrixIsConst = time ~= 0.0 and dtChanged == false
	assTuner:set_matrix_is_const(matrixIsConst)
	AssembleLinearOperatorRhsAndSolution(linOp, u, b)
	
	-- apply linear solver
	ilu:set_disable_preprocessing(matrixIsConst)
	if ApplyLinearSolver(linOp, u, b, solver) == false then
		print("Could not apply linear solver.")
		exit()
	end
	
	-- update to new time
	time = solTimeSeries:time(0) + curr_dt
	
	-- vtk output
	if generateVTKoutput then
		if math.abs(time/pstep - math.floor(time/pstep+0.5)) < 1e-5 then 
			out:print(outDir.."vtk/solution", u, math.floor(time/pstep+0.5), time)
		end
	end
	
	take_measurement(u, time, "syn", measFcts, outDir.."meas/sol")
	
	-- updte time series (reuse memory)
	oldestSol = solTimeSeries:oldest()
	VecScaleAssign(oldestSol, 1.0, u)
	solTimeSeries:push_discard_oldest(oldestSol, time)
	
	-- increment check-back counter
	cb_counter[lv] = cb_counter[lv] + 1

	print("++++++ POINT IN TIME " .. math.floor((time)/curr_dt+0.5)*curr_dt .. "  END ++++++")
end

-- end timeseries, produce gathering file
if generateVTKoutput then 
	out:write_time_pvd(outDir.."vtk/solution", u) 
end

currentOutFile:close()

	
