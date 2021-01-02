--------------------------------------------------------------------------------
-- Script for PNP simulation on reconstructed spine geometry                  --
-- (with 1D extensions and ER).                                               --
--                                                                            --
-- Author: Markus Breit                                                       --
-- Date:   2019-01-21                                                         --
--------------------------------------------------------------------------------

-- load pre-implemented lua functions
ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")

-- init with dimension and algebra
InitUG(3, AlgebraType("CPU", 5))


-----------------------------------------------------
--  parameters steering the simulation and output  --
-----------------------------------------------------
-- choice of grid
gridName = util.GetParam("-grid", "grids/spineER.ugx")

-- total refinements
numChargeRefs = util.GetParamNumber("-numChargeRefs", 0)
numGlobRefs = util.GetParamNumber("-numGlobRefs", 0)

-- time stepping
startTime = util.GetParamNumber("-start", 0.0, "start time")
endTime = util.GetParamNumber("-end", 0.02, "end time")
dt = util.GetParamNumber("-dt", 1e-4, "time step size")
dtmin = util.GetParamNumber("-dtmin", 1e-7, "minimal admissible time step size")
dtmax = util.GetParamNumber("-dtmax", 1e-2, "maximal admissible time step size")

-- choose outfile directory
outDir = util.GetParam("-outName", "PNP_full")
outDir = outDir.."/"

-- error tolerance and number of stages for Limex iteration
toleratedError = util.GetParamNumber("-tol", 0.001)
nstages = util.GetParamNumber("-nst", 2)

-- additional verbosity, e.g. linear solver output
verbose	= util.HasParamOption("-verbose")

-- allowed processor imbalance
imbFactor = util.GetParamNumber("-imb", 1.05, "imbalance factor")

-- specify -vtk to generate vtk output
generateVTKoutput = util.HasParamOption("-vtk")
pstep = util.GetParamNumber("-pstep", dt, "plotting interval")

-- which upwind method to use
upwindMethod = util.GetParam("-upwind", "FULL")
upwindMethod = string.upper(upwindMethod)
if upwindMethod ~= "PNP" and upwindMethod ~= "FULL" and upwindMethod ~= "NO" then
	print("Not a viable upwind method: '" .. upwindMethod .. "'.")
	exit()
end

-- instationary calculation?
bInstationary = util.HasParamOption("-instat")

-- profiling?
doProfiling = util.HasParamOption("-profile")
SetOutputProfileStats(doProfiling)

-- debugging gmg?
gmgDebug = util.HasParamOption("-gmgDebug")

-- check whether to import solution
bImportSol = util.HasParamOption("-import")
chkptFile = util.GetParam("-chkptFile", "noChkptFileProvided") -- without .lua !!


-------------------------
--  problem constants  --
-------------------------
Phi				=	"Phi"
ion_species		=	{   "K"   ,  "Cl"  ,  "Na"  ,   "A" }
valency			=	{     1   ,   -1   ,    1   ,   -1  }
conc_out      	= 	{     4   ,  123   ,  145   ,   26  }   -- in mol/m^3 == mM
conc_in			= 	{   155   ,    5   ,   12   ,  162  }   -- in mol/m^3 == mM
diff_const		=	{  1.96e-9, 2.03e-9, 1.33e-9, 2.0e-9}   -- in m^2/s
reversal_pot	=	{ -0.07   , -0.07  , -0.07	, -0.07	}	-- in V
spec_cond		=	{   0.5   ,    0   ,  0.5   ,    0  }	-- in S/m^2
PSD_channel_dens=   {  1e-4   ,    0   , 1e-4   ,    0  }	-- in 1

F = 96485  -- in C / mol
R = 8.31451  -- in J / (mol K)
T = 298.15  -- in K
vac_perm = 8.8541878e-12  -- in C / (Vm)

rel_perm_mem = 4.0
rel_perm_rest = 80.0

pot_out = 0.0
pot_in = -0.070

surfCh = util.GetParamNumber("-surfCh", 0.0)

--geometry
radius = 4e-7 -- 400nm
mem_thickness = radius*math.log(1.0 + 1e-8/radius)

-- membrane capacitances
total_conc = 0
for i = 1, 4 do
	total_conc = total_conc + valency[i]*valency[i]*conc_in[i]
end
cap_factor = rel_perm_mem * vac_perm / mem_thickness / total_conc
spec_cap = {}  -- in F/m^2
for i = 1, 4 do
	spec_cap[i] = valency[i]*valency[i]*conc_in[i] * cap_factor
end


-- additional scaling of equations
pot_fac = 1.0/F

-- dimless parameters --
conc_ref = 100.0
pot_ref = R*T/F
diff_ref = 2e-9
len_ref = 1e-6
time_ref = len_ref*len_ref/diff_ref
cond_ref = conc_ref * diff_ref / (pot_ref * len_ref) * F
cap_ref = conc_ref * len_ref / pot_ref * F
flux_ref = diff_ref * conc_ref / len_ref

dimless_startTime = startTime / time_ref
dimless_endTime = endTime / time_ref
dimless_dt = dt / time_ref
dimless_dtmin = dtmin / time_ref
dimless_dtmax = dtmax / time_ref

dimless_conc_out = {}
dimless_conc_in = {}
dimless_diff_const = {}
dimless_reversal_pot = {}
dimless_spec_cond = {}
dimless_spec_cap = {}

for i = 1,4 do
	dimless_conc_out[i] = conc_out[i] / conc_ref
	dimless_conc_in[i] = conc_in[i] / conc_ref
	dimless_diff_const[i] = diff_const[i] / diff_ref
	dimless_reversal_pot[i] = reversal_pot[i] / pot_ref
	dimless_spec_cond[i] = spec_cond[i] / cond_ref
	dimless_spec_cap[i] = spec_cap[i] / cap_ref
end

dimless_pot_out = pot_out / pot_ref
dimless_pot_in = pot_in / pot_ref

dimless_surfCh = surfCh / (conc_ref*len_ref) * pot_fac

dimless_perm_mem = rel_perm_mem * vac_perm * pot_ref / (len_ref*len_ref*conc_ref) * pot_fac
dimless_perm_rest = rel_perm_rest * vac_perm * pot_ref / (len_ref*len_ref*conc_ref) * pot_fac

dimless_radius = radius / len_ref
dimless_mem_thickness = mem_thickness / len_ref

dim_scales = {conc_ref, conc_ref, conc_ref, conc_ref, pot_ref}

print("Reference time        : " .. time_ref)
print("Reference length      : " .. len_ref)
print("Reference potential   : " .. pot_ref)
print("Reference conentration: " .. conc_ref)


----------------------------
--  setup discretization  --
----------------------------
dom = Domain()
dom:create_additional_subset_handler("projSH")
LoadDomain(dom, gridName)

neededSubsets = {"in", "mem", "psd", "er", "mem_in", "psd_in", "mem_out", "psd_out",
	"mem_er", "bnd_mem", "intfLeft_constr", "intfRight_constr",
	"intfLeft_node1D", "intfRight_node1D", "intfLeft_nodehD", "intfRight_nodehD",
	"extLeft", "extRight", "bnd_extLeft", "bnd_extRight", "useless"}
ug_assert(util.CheckSubsets(dom, neededSubsets), "Something wrong with required subsets.")

-- create approximation space
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct(ion_species, "Lagrange", 1)	-- add species
approxSpace:add_fct(Phi, "Lagrange", 1)			-- add electric potential
approxSpace:init_top_surface()
approxSpace:init_levels()


-- save unknowns as string and vector
all_unknowns = {}
all_unknowns_string = ""
all_ions_string = ""
for i=1, #ion_species do
	all_unknowns[i] = ion_species[i]
	all_unknowns_string = all_unknowns_string .. ion_species[i]..","
	all_ions_string = all_ions_string .. ", " .. ion_species[i]
end
all_unknowns[#all_unknowns + 1] = Phi
all_unknowns_string = all_unknowns_string .. Phi
all_ions_string = string.sub(all_ions_string, 3)

-- name some subdomains
local intra = "in, extLeft, extRight"
local filMem = "er, mem, psd, bnd_mem"
local chSurfIn = "mem_in, psd_in"
local intfcs = "intfLeft_nodehD, intfLeft_node1D, intfLeft_constr, " ..
               "intfRight_nodehD, intfRight_node1D, intfRight_constr"

-- create domain disc
domainDisc = DomainDiscretization(approxSpace)

-- upwind
upwind = PNPUpwind()
upwind:set_alpha(1e-3)
if upwindMethod == "NO" then
	upwind = NoUpwind()
elseif upwindMethod == "FULL" then
	upwind = FullUpwind()
end


-- inside: fully coupled Nernst-Planck equations
sourceCharge = ScaleAddLinkerNumber()
potentialDisc = ConvectionDiffusionFV1(Phi, "in")
potentialDisc:set_stationary()
potentialDisc:set_diffusion(dimless_perm_rest)
potentialDisc:set_source(sourceCharge)
domainDisc:add(potentialDisc)
	
for i = 1, #ion_species do
	-- compute velocity
	local vel = ScaleAddLinkerVector()
	vel:add(-valency[i]*F/(R*T)*pot_ref*dimless_diff_const[i], potentialDisc:gradient())

	-- equation for species
	local speciesDisc = ConvectionDiffusionFV1(ion_species[i], "in")
	speciesDisc:set_diffusion(dimless_diff_const[i])
	speciesDisc:set_velocity(vel)
	speciesDisc:set_upwind(upwind)
	domainDisc:add(speciesDisc)
	
	-- contribution to charge source
	sourceCharge:add(valency[i]*F*pot_fac, speciesDisc:value())			
end

-- membrane, ER: only electro-static problem for potential, ohmic conductivity on Ranvier nodes
estatPotentialDisc = ConvectionDiffusionFV1(Phi, "mem, psd, er")
estatPotentialDisc:set_stationary()
estatPotentialDisc:set_diffusion(dimless_perm_mem)
domainDisc:add(estatPotentialDisc)

memSpec = DirichletBoundary()
for i = 1, #ion_species do
	memSpec:add(0.0, ion_species[i], "er, mem, psd, bnd_mem")
end
--memSpec:add(dimless_conc_in[2], ion_species[2], "psd")
--memSpec:add(dimless_conc_in[4], ion_species[4], "psd")
domainDisc:add(memSpec)


-- add surface charges on filaments and membrane
if dimless_surfCh ~= 0 then
	local NeumannBND = UserFluxBoundaryFV1(Phi, chSurfIn)
	NeumannBND:set_flux_function(dimless_surfCh)
	NeumannBND:set_stationary()
	domainDisc:add(NeumannBND)
end


-- PNP 1D elem disc
dimFacVol = math.pi * dimless_radius * dimless_radius
dimFacSurf = 2.0 * math.pi * dimless_radius
	
sourceCharge1d = ScaleAddLinkerNumber()
potentialDisc1d = ConvectionDiffusionFV1(Phi, "extLeft, extRight")
potentialDisc1d:set_stationary()
potentialDisc1d:set_diffusion(dimFacVol*dimless_perm_rest)
potentialDisc1d:set_source(dimFacVol*sourceCharge1d)
domainDisc:add(potentialDisc1d)

for i = 1, #ion_species do
	-- compute "velocity"
	local vel = ScaleAddLinkerVector()
	vel:add(-dimFacVol*dimless_diff_const[i]*valency[i]*F/(R*T)*pot_ref, potentialDisc1d:gradient())

	-- mass
	local mass = ScaleAddLinkerNumber()
	mass:add(dimFacSurf*dimless_spec_cap[i]/valency[i], potentialDisc1d:value())

	-- source
	local source = ScaleAddLinkerNumber()
	source:add(-dimFacSurf*dimless_spec_cond[i]/valency[i], potentialDisc1d:value())
	source:add(dimFacSurf*dimless_spec_cond[i]/valency[i], dimless_reversal_pot[i])

	-- equation for species
	local speciesDisc = ConvectionDiffusionFV1(ion_species[i], "extLeft, extRight")
	speciesDisc:set_upwind(upwind)
	speciesDisc:set_diffusion(dimFacVol*dimless_diff_const[i])
	speciesDisc:set_velocity(vel)
	speciesDisc:set_mass(mass)
	speciesDisc:set_mass_scale(dimFacVol)
	if spec_cond[i] ~= 0 then
		speciesDisc:set_source(source)
	end
	
	domainDisc:add(speciesDisc)
	
	-- contribution to charge source
	sourceCharge1d:add(valency[i]*F*pot_fac, speciesDisc:value())			
end


-- hD/1D interfaces
PotIntfLeft = AdditiveInterface1D(Phi, "intfLeft_constr", "intfLeft_nodehD", "intfLeft_node1D", {-1, 0, 0})
PotIntfRight = AdditiveInterface1D(Phi, "intfRight_constr", "intfRight_nodehD", "intfRight_node1D", {1, 0, 0})
IonIntfLeft = MultiplicativeInterface1D(all_ions_string, "intfLeft_constr", "intfLeft_nodehD", "intfLeft_node1D", {-1, 0, 0})
IonIntfRight = MultiplicativeInterface1D(all_ions_string, "intfRight_constr", "intfRight_nodehD", "intfRight_node1D", {1, 0, 0})
domainDisc:add(PotIntfLeft)
domainDisc:add(PotIntfRight)
domainDisc:add(IonIntfLeft)
domainDisc:add(IonIntfRight)


-- Dirichlet boundaries
diriBnds = DirichletBoundary()
diriBnds:add(dimless_pot_out, Phi, "mem_out, psd_out, useless")
for i = 1, #ion_species do
	diriBnds:add(dimless_conc_out[i], ion_species[i], "mem_out, psd_out, useless")
end
domainDisc:add(diriBnds)

-- extra Dirichlet bnds for the stationary case
statDiriBnds = DirichletBoundary()
statDiriBnds:add(dimless_pot_in, Phi, "bnd_extLeft, bnd_extRight")
for i = 1, #ion_species do
	statDiriBnds:add(dimless_conc_in[i], ion_species[i], "bnd_extLeft, bnd_extRight")
end
--statDiriBnds:add(dimless_conc_in[1], ion_species[1], "psd")
--statDiriBnds:add(dimless_conc_in[3], ion_species[3], "psd")
domainDisc:add(statDiriBnds)

-- constraints for adaptivity
hangingConstraint = SymP1Constraints()
domainDisc:add(hangingConstraint)


-- create stationary operator from domain discretization
op_stat = AssembledOperator(domainDisc)


-- time discretization
timeDisc = ThetaTimeStep(domainDisc)
timeDisc:set_theta(1.0) -- 1.0 is implicit Euler


------------------
-- solver setup --
------------------
-- debug writer
dbgWriter = GridFunctionDebugWriter(approxSpace)
dbgWriter:set_base_dir(outDir)
dbgWriter:set_vtk_output(false)

-- linear convergence check
convCheck = ConvCheck()
convCheck:set_minimum_defect(1e-13)
convCheck:set_reduction(1e-06)
convCheck:set_verbose(verbose)
--[[
convCheck = CompositeConvCheck(approxSpace, 10000, 1e-50, 1e-04)
convCheck:set_component_check(Phi, 1e-13, 1e-06)
convCheck:set_component_check(ion_species, 1e-13, 1e-06)
--convCheck:set_rest_check(1e-13, 1e-04)
convCheck:set_verbose(verbose)
convCheck:set_adaptive(true)
--]]

linearSolver = BiCGStab()--LinearSolver()
linearSolver:set_convergence_check(convCheck)
if gmgDebug then
	linearSolver:set_debug(dbgWriter)
end

gmg = GeometricMultiGrid(approxSpace)
gmg:set_discretization(domainDisc)
gmg:set_base_level(0)
gmg:set_gathered_base_solver_if_ambiguous(true)
gmg:set_base_solver(SuperLU())

smoother = GaussSeidel()
smoother:enable_consistent_interfaces(true)
--smoother:set_sort(true)

gmg:set_smoother(smoother)
gmg:set_smooth_on_surface_rim(true)
gmg:set_cycle_type(1)
gmg:set_num_presmooth(3)
gmg:set_num_postsmooth(3)
--gmg:set_rap(true) -- causes error in base solver!!

if gmgDebug then
	gmg:set_debug(dbgWriter)
end

linearSolver:set_preconditioner(gmg)
convCheck:set_maximum_steps(100)


--- non-linear solver ---
newtonConvCheck = CompositeConvCheck(approxSpace, 10, 1e-22, 1e-10)
 -- minDefect values painstakingly estimated: do not change!
newtonConvCheck:set_component_check(Phi, 1e-9, 1e-10)
newtonConvCheck:set_component_check(ion_species, 2e-12*math.pow(1, 2), 1e-10)
newtonConvCheck:set_verbose(true)
newtonConvCheck:set_time_measurement(true)
newtonConvCheck:set_adaptive(true)

-- Newton solver
newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(linearSolver)
newtonSolver:set_convergence_check(newtonConvCheck)
--newtonSolver:set_debug(dbgWriter)


---------------------------------
--  refinement & distribution  --
---------------------------------
-- in parallel environments: use a load balancer to distribute the grid
-- actual refinement and load balancing after setup of disc.
balancer.partitioner = "parmetis"
balancer.staticProcHierarchy = true
balancer.firstDistLvl   = 0
balancer.redistSteps    = 2
balancer.maxDistLvl     = 4
balancer.firstDistProcs = 48
balancer.redistProcs    = 64
balancer.parallelElementThreshold = 8

balancer.imbalanceFactor = imbFactor
balancer.ParseParameters()
balancer.PrintParameters()

loadBalancer = balancer.CreateLoadBalancer(dom)


-- refinement and distribution --
refiner = HangingNodeDomainRefiner(dom)

-- configure distro adjuster
ida = PNPDistroManager(approxSpace)
ida:add_interface(PotIntfLeft)
ida:add_interface(PotIntfRight)
ida:add_interface(IonIntfLeft)
ida:add_interface(IonIntfRight)

ida:set_threshold_ratio(0.25) -- anisotropy threshold (elements below will be glued to their neighbors)

cdgm = ClusteredDualGraphManager()
cdgm:add_unificator(SiblingUnificator())
cdgm:add_unificator(ida)

-- add distro adjuster to distributed grid manager to be considered when distributing
set_distro_adjuster(dom, ida)
if loadBalancer ~= nil and balancer.partitioner == "parmetis" then
	balancer.defaultPartitioner:set_dual_graph_manager(cdgm)
end

-- set interface refinement adjuster
irma = InterfaceRefMarkAdjuster()
irma:set_subset_handler(dom:subset_handler())
irma:add_interfaces({PotIntfLeft, PotIntfRight, IonIntfLeft, IonIntfRight})
add_interface_ref_mark_adjuster(refiner, irma)

-- set extension refinement adjuster
erma = ExtensionRefMarkAdjuster(dom, {1,0,0}, "useless")
add_extension_ref_mark_adjuster(refiner, erma)

-- initial rebalancing
if loadBalancer ~= nil then
	balancer.qualityRecordName = "init"
	balancer.Rebalance(dom, loadBalancer)
end


-- charge refinements
strat = SurfaceMarking(dom)
--strat:add_surface("mem_in", "in")
strat:add_surface("psd_in", "in")

for i = 1, numChargeRefs do
	strat:mark_without_error(refiner, approxSpace)
	refiner:refine()

	-- rebalance
	if loadBalancer ~= nil then
		balancer.qualityRecordName = "chrgRef " .. i
		balancer.Rebalance(dom, loadBalancer)
	end
end

-- global refinements	
for i = 1, numGlobRefs do
	mark_global(refiner, dom)
	refiner:refine()
	
	-- rebalance
	if loadBalancer ~= nil then
		balancer.qualityRecordName = "globRef " .. i
		balancer.Rebalance(dom, loadBalancer)
	end
end

--TestDomainInterfaces(dom)

print()
approxSpace:init_levels()	-- re-init for new levels
approxSpace:init_surfaces()
approxSpace:init_top_surface()
approxSpace:print_statistic()
print()
print(dom:domain_info():to_string())

if loadBalancer ~= nil then
	loadBalancer:estimate_distribution_quality()
	loadBalancer:print_quality_records()
end

-- in highly refined case, this takes forever AND has a REALLY big result
--SaveGridHierarchyTransformed(dom:grid(), dom:subset_handler(), outDir.."grid/refined_grid_hierarchy_lv".. numGlobRefs+numChargeRefs .."_p" .. ProcRank() .. ".ugx", 5.0)
--SaveParallelGridLayout(dom:grid(), outDir.."grid/parallel_grid_layout_lv".. numGlobRefs+numChargeRefs .."_p"..ProcRank()..".ugx", 5.0)

-- we can only rescale AFTER all refinements (otherwise, projectors are scaled incorrectly)
if len_ref ~= 1e-6 then	 -- domain is saved in units of um
	RescaleDomain(dom, 1e-6/len_ref)
end


-- we can only sort the dofs AFTER all refinements (otherwise, the sorting is not performed on highest levels)
print("\nordering using Cuthill-McKee\n")
OrderCuthillMcKee(approxSpace, true) -- this is essential (at least in one 2d spine problem)!

-- update interfaces
PotIntfLeft:update()
PotIntfRight:update()
IonIntfLeft:update()
IonIntfRight:update()


----------------------
--  Initial values  --
----------------------
u = GridFunction(approxSpace)
u:set(0.0)
u_dispose = u:clone()


InterpolateInner(dimless_pot_in, u, "Phi", intra .. ", er, mem_er, " .. intfcs .. ", bnd_extLeft, bnd_extRight")

d_phi = 0.0
if surfCh ~= 0.0 then d_phi = 0.5*math.abs(dimless_pot_in) end
InterpolateInner(dimless_pot_in-d_phi, u, "Phi", chSurfIn..", mem, psd, bnd_mem")
 
for i = 1, #ion_species do
	InterpolateInner(dimless_conc_out[i], u, ion_species[i], "mem_out, psd_out")
	InterpolateInner(dimless_conc_in[i] , u, ion_species[i], intra .. ", mem_er, "
					 .. intfcs .. ", bnd_extLeft, bnd_extRight")
end

cations = {[1] = ion_species[1], [3] = ion_species[3]}
anions = {[2] = ion_species[2], [4] = ion_species[4]}

factor = 1.0
if surfCh ~= 0.0 then factor = 2.0 end
for i,cat in pairs(cations) do
	InterpolateInner(factor*dimless_conc_in[i], u, cat, chSurfIn)
end
if surfCh ~= 0.0 then factor = 0.5 end
for i,an in pairs(anions) do
	InterpolateInner(factor*dimless_conc_in[i], u, an, chSurfIn)
end


-------------------------------
--  vtk output preparations  --
-------------------------------
out = VTKOutput()
out:select_all(true)
out_flux = VTKOutput()

solAdjuster = Domain1dSolutionAdjuster()
solAdjuster:add_constrained_subset("useless")
solAdjuster:add_constrainer_subset("extLeft")
solAdjuster:add_constrainer_subset("intfLeft_node1D")
solAdjuster:add_constrainer_subset("extRight")
solAdjuster:add_constrainer_subset("intfRight_node1D")
solAdjuster:set_sorting_direction({1,0,0})

-- init neck recorder
nr = NeckRecorder(approxSpace)
nr:set_cytosolic_subset("in")
nr:add_measurement_zone(0.5*1e-6/len_ref, "base")
nr:add_measurement_zone(1.2*1e-6/len_ref, "mid")
nr:add_measurement_zone(1.9*1e-6/len_ref, "head")
nr:set_diffusion_constants({dimless_diff_const[1], dimless_diff_const[3],
							dimless_diff_const[2], dimless_diff_const[4]})
nr:set_convection_constants(
	{-dimless_diff_const[1]*valency[1]*F/(R*T)*pot_ref,
     -dimless_diff_const[3]*valency[3]*F/(R*T)*pot_ref,
     -dimless_diff_const[2]*valency[2]*F/(R*T)*pot_ref,
     -dimless_diff_const[4]*valency[4]*F/(R*T)*pot_ref})
nr:set_record_individual_currents(true)
nr:set_temperature(T)
nr:set_upwind(upwind)


-------------------------------------
--  calculate stationary solution  --
-------------------------------------
if bImportSol then
	cp = util.ReadCheckpoint(u, chkptFile)
else
	newtonSolver:init(op_stat)
	if not newtonSolver:apply(u) then
		print("Newton Solver failed.")
		if generateVTKoutput then
			ScaleGF(u_dispose, u, dim_scales)
	  	 	solAdjuster:adjust_solution(u_dispose)
			out:print(outDir.."vtk/failedSolution", u_dispose, 0, 0.0)
		end
		if doProfiling then
			WriteProfileData(outDir .."pd.pdxml")
		end
		exit()
	end

	-- export stationary solution
	--util.WriteCheckpoint(u, 0, {dimless_time = 0.0, step = 0}, outDir .. "chkpt")
end


-- write stationary solution
ScaleGF(u_dispose, u, dim_scales)
if generateVTKoutput then
	solAdjuster:adjust_solution(u_dispose)
	out:print(outDir.."vtk/solution", u_dispose, 0, 0.0)
end
--take_measurement(u_dispose, 0.0, "psd_in", "Phi, K, Na, Cl, A", outDir.."meas/sol")

-- record initial values
nr:record_current(outDir .."meas/neck_current.dat", 0.0, u, F * diff_ref*conc_ref * len_ref)
nr:record_potential(outDir .."meas/neck_pot.dat", 0.0, u, pot_ref)
nr:record_concentrations(outDir .."meas/neck_concs.dat", 0.0, u, conc_ref)


-------------------------
--  compute time loop  --
-------------------------
-- continue only if -instat flag is set
if bInstationary == true then


-- synaptic activation --
local t_on = 0.0
local tau = 2.5e-3  -- in s
local t_on_dimless = t_on / time_ref
local tau_dimless = tau / time_ref
local period_dimless = 0.01 / time_ref
local tau_dimless_nmdar = 0.15 / time_ref

function channelOpen(t)
	-- strong NMDAR
	--return math.exp(-t/tau_dimless_nmdar)
	
	-- tetanic stimulation
	--local tmp = math.sin(math.pi*t/period_dimless)
	--return tmp*tmp
	
	-- AMPAR
	return (t - t_on_dimless) / tau_dimless * math.exp(-(t - t_on_dimless - tau_dimless) / tau_dimless)
end

---[[
vmAtSyn = pot_in
kInAtSyn = dimless_conc_in[1]
naInAtSyn = dimless_conc_in[3]

psdInArea = Integral(1.0, u, "psd_in")
psdOutArea = Integral(1.0, u, "psd_out")

function psd_flux_k_in(x,y,z,t,si)
	local exp = math.exp(F/(R*T) * vmAtSyn)
	return - channelOpen(t) * PSD_channel_dens[1] * dimless_diff_const[1] / 1e-8 * len_ref *
		F/(R*T) * vmAtSyn * (dimless_conc_out[1] - kInAtSyn * exp) / (1.0 - exp)
end

function psd_flux_na_in(x,y,z,t,si)
	local exp = math.exp(F/(R*T) * vmAtSyn)
	return - channelOpen(t) * PSD_channel_dens[3] * dimless_diff_const[3] / 1e-8 * len_ref *
		F/(R*T) * vmAtSyn * (dimless_conc_out[3] - naInAtSyn * exp) / (1.0 - exp)
end

PSD_NeumannKInner = UserFluxBoundaryFV1(ion_species[1], "psd_in")
PSD_NeumannKInner:set_flux_function("psd_flux_k_in")
domainDisc:add(PSD_NeumannKInner)

PSD_NeumannNaInner = UserFluxBoundaryFV1(ion_species[3], "psd_in")
PSD_NeumannNaInner:set_flux_function("psd_flux_na_in")
domainDisc:add(PSD_NeumannNaInner)
--]]

--[[

lunChannelOpen = LuaUserNumber("channelOpen")

function channelOpenTensor(x,y,z,t,si)
	local p = channelOpen(x,y,z,t,si)
	return p, 0, 0,
	       0, p, 0,
	       0, 0, p
end
lumChannelOpen = LuaUserMatrix("channelOpenTensor")

psdInArea = Integral(1.0, u, "psd_in")
areaFactor = 0.09275 / psdInArea

psdVel1 = ScaleAddLinkerVector()
psdVel1:add(-valency[1]*F/(R*T)*pot_ref*dimless_diff_const[1]*PSD_channel_dens[1]*areaFactor*lunChannelOpen, estatPotentialDisc:gradient())
psdSpeciesDisc1 = ConvectionDiffusionFV1(ion_species[1], "psd")
psdSpeciesDisc1:set_diffusion(dimless_diff_const[1]*PSD_channel_dens[1]*areaFactor*lumChannelOpen)
psdSpeciesDisc1:set_velocity(psdVel1)
psdSpeciesDisc1:set_upwind(upwind)
domainDisc:add(psdSpeciesDisc1)

psdVel3 = ScaleAddLinkerVector()
psdVel3:add(-valency[3]*F/(R*T)*pot_ref*dimless_diff_const[3]*PSD_channel_dens[3]*areaFactor*lunChannelOpen, estatPotentialDisc:gradient())
psdSpeciesDisc3 = ConvectionDiffusionFV1(ion_species[3], "psd")
psdSpeciesDisc3:set_diffusion(PSD_channel_dens[3]*PSD_channel_dens[3]*areaFactor*lumChannelOpen)
psdSpeciesDisc3:set_velocity(psdVel3)
psdSpeciesDisc3:set_upwind(upwind)
domainDisc:add(psdSpeciesDisc3)
--]]

-- remove Dirichlet bnd from extensions
-- and replace by Neumann-0 (which happens implicitly by doing nothing)
-- (GMG transfer operators must be updated!)
domainDisc:remove(statDiriBnds)
gmg:force_reinit()


convCheck:set_reduction(1e-6)

-- create instationary operator from time discretization
op = AssembledOperator(timeDisc)


if bImportSol then
	dimless_time = cp.myData.dimless_time
	step = cp.myData.step
else
	dimless_time = dimless_startTime
	step = 0
end




--  LIMEX setup --
-- create instationary solver
newtonSolver = LimexNewtonSolver()
newtonSolver:set_linear_solver(linearSolver)
newtonSolver:init(op)

stageNSteps = {}    -- number of time steps for each stage
for i = 1, nstages do stageNSteps[i] = i end

limex = LimexTimeIntegrator(nstages)
for i = 1, nstages do
	limex:add_stage(stageNSteps[i], newtonSolver, domainDisc)
end

limex:set_tolerance(toleratedError)
limex:set_time_step(dimless_dt)
limex:set_dt_min(dimless_dtmin)
limex:set_dt_max(dimless_dtmax)
limex:set_increase_factor(2.0)
limex:set_reduction_factor(0.1)
limex:set_stepsize_greedy_order_factor(1)
limex:set_stepsize_safety_factor(0.25)

-- GridFunction error estimator (relative norm)
errorEvalPhi = H1ComponentSpace("Phi")--H1SemiComponentSpace("Phi", 3, dimless_perm_rest)  -- function name, order, weight
errorEvalK = H1ComponentSpace("K")
errorEvalNa = H1ComponentSpace("Na")
errorEvalCl = H1ComponentSpace("Cl")
errorEvalA = H1ComponentSpace("A")

limexEstimator = GridFunctionEstimator()--toleratedError)
limexEstimator:add(errorEvalPhi)
limexEstimator:add(errorEvalK)
limexEstimator:add(errorEvalNa)
limexEstimator:add(errorEvalCl)
limexEstimator:add(errorEvalA)
limex:add_error_estimator(limexEstimator)


-- for post-processing after each time step
function postTimestepActions(step, time, dt)
	local curSol = preAndPostActions:get_current_solution()

	-- plot solution every pstep seconds
	ScaleGF(u_dispose, curSol, dim_scales)  -- scale to usual units
	if generateVTKoutput then
		solAdjuster:adjust_solution(u_dispose)
		out:print(outDir.."vtk/solution", u_dispose, step, time*time_ref)
	end
	--take_measurement(u_dispose, time*time_ref, "psd_in", "Phi, K, Na, Cl, A", outDir.."meas/sol")

	-- record current and potential
	nr:record_current(outDir .."meas/neck_current.dat", time*time_ref, curSol, F * diff_ref*conc_ref * len_ref)
	nr:record_potential(outDir .."meas/neck_pot.dat", time*time_ref, curSol, pot_ref)
	nr:record_concentrations(outDir .."meas/neck_concs.dat", time*time_ref, curSol, conc_ref)
	
	-- export stationary solution
	--util.WriteCheckpoint(curSol, time*time_ref, {dimless_time = time, step = step}, outDir .. "chkpt")
	
	print("Current (real) time: " .. time*time_ref .. ",   last dt: " .. dt*time_ref)
	
	return 0
end

function preTimestepActions(step, time, dt)
	local curSol = preAndPostActions:get_current_solution()
	vmAtSyn = pot_ref * (Integral(curSol, Phi, "psd_in") / psdInArea)
	
	---[[
	kInAtSyn = Integral(curSol, "K", "psd_in") / psdInArea
	naInAtSyn = Integral(curSol, "Na", "psd_in") / psdInArea
	--]]
	
	return 0
end

preAndPostActions = LuaCallbackObserver()
preAndPostActions:set_callback_pre("preTimestepActions")
preAndPostActions:set_callback_post("postTimestepActions")
limex:attach_observer(preAndPostActions)


-- solve problem
limex:apply(u, dimless_endTime, u, dimless_time)



end -- if bInstationary


-- end timeseries, produce gathering file
if generateVTKoutput then
	out:write_time_pvd(outDir .. "vtk/solution", u)
end


if doProfiling then
	WriteProfileData(outDir .."pd.pdxml")
end


