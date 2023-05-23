"""
21600 points
102448 bonds

Note that recent changes on Peridynamics@main haven't been tagged
"""

using KitPD
using ProgressMeter

lx = ly = 0.05
lz = 0.005
ny = 60
Δx = ly / ny
pc = PD.PointCloud(lx, ly, lz, Δx)

mat = PD.BondBasedMaterial(;
	horizon = 3.015Δx,
	rho = 7850.0,
	E = 210e9,
	Gc = 1000.0,
)

begin
	cracklength = 0.5 * lx
	precrack_set_a = findall(
		(pc.position[2, :] .>= 0) .&
		(pc.position[2, :] .< 12 * Δx) .& # y ∈ [0, 0.01]
		(pc.position[1, :] .<= -lx / 2 + cracklength), # x < 0
	)
	precrack_set_b = findall(
		(pc.position[2, :] .<= 0) .&
		(pc.position[2, :] .> -12 * Δx) .& # y ∈ [-0.01, 0]
		(pc.position[1, :] .<= -lx / 2 + cracklength),
	)
	precracks = [PD.PreCrack(precrack_set_a, precrack_set_b)]
end

begin
	bc_set_top = findall(pc.position[2, :] .> ly / 2 - 5.1 * Δx)
	bc_set_bottom = findall(pc.position[2, :] .< -ly / 2 + 5.1 * Δx)
	bc_top = PD.VelocityBC(t -> 0.1, bc_set_top, 2)
	bc_bottom = PD.VelocityBC(t -> -0.1, bc_set_bottom, 2)
	boundary_conditions = [bc_top, bc_bottom]
end

td = PD.TimeDiscretization(2000)

begin
	simulation_name = "Test"
	resfolder = joinpath(@__DIR__, "results", simulation_name)
	mkpath(resfolder)
	es = PD.ExportSettings(resfolder, 100)
end

sim = PD.PDSingleBodyAnalysis(;
	name = simulation_name,
	pc = pc,
	mat = mat,
	precracks = precracks,
	bcs = boundary_conditions,
	td = td,
	es = es,
)

body = PD.create_simmodel(sim.mat, sim.pc)
for precrack in sim.precracks
	PD.define_precrack!(body, precrack)
end

### test
findall(body.bond_failure .== 0)
body.bond_data[7874]
163 in precrack_set_b && 181 in precrack_set_a

findall(0 .< body.n_active_family_members[:, 1] .< body.n_family_members)
###

PD.update_thread_cache!(body)
PD.calc_damage!(body)

### test
findall(body.damage .== 1)
###

if sim.td.Δt < 0.0 && sim.td.alg !== :dynrelax
	sim.td.Δt = PD.calc_stable_timestep(body, sim.mat.rho, sim.mat.K, sim.mat.δ)
end

PD.export_vtk(body, sim.es.resultfile_prefix, 0, 0.0)

p = Progress(sim.td.n_timesteps;
	dt = 1,
	desc = "Time integration... ",
	barlen = 30,
	color = :normal,
	enabled = !PD.is_logging(stderr),
)
Δt½ = 0.5 * sim.td.Δt
begin
    for t in 1:2000#sim.td.n_timesteps
        time = t * sim.td.Δt
        PD.update_velhalf!(body, Δt½)
        PD.apply_bcs!(body, sim.bcs, time)
        PD.update_disp_and_position!(body, sim.td.Δt)
        PD.compute_forcedensity!(body, sim.mat)
        PD.update_thread_cache!(body)
        PD.calc_damage!(body)
        PD.compute_equation_of_motion!(body, Δt½, sim.mat.rho)
        if mod(t, sim.es.exportfreq) == 0
            PD.export_vtk(body, sim.es.resultfile_prefix, t, time)
        end
        next!(p)
    end
finish!(p)
end
