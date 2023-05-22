"""
21600 points
102448 bonds
"""

using KitPD

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






