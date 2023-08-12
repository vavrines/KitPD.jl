"""
Updated crack example based on Peridynamics.jl v0.2.0
Differences include (not limited to):
- ContinuumBasedMaterial
- VelocityVerlet
- init_time_discretization!
- compute_equation_of_motion!
"""

using KitPD
using KitBase.ProgressMeter

cd(@__DIR__)

l = 1.0
Δx = l / 50
pc = PD.PointCloud(l, l, 0.1l, Δx)

δ = 3.015 * Δx
mat = PD.ContinuumBasedMaterial(horizon = δ, rho = 8e-6, E = 2.1e5, nu = 0.25, Gc = 2.7)

a = 0.5 * l
set_a = findall(p -> p[1] ≤ -l / 2 + a && 0 ≤ p[2] ≤ 2δ, eachcol(pc.position))
set_b = findall(p -> p[1] ≤ -l / 2 + a && -2δ ≤ p[2] < 0, eachcol(pc.position))
precrack = PD.PreCrack(set_a, set_b)

set_top = findall(p -> p[2] > l / 2 - Δx, eachcol(pc.position))
set_bottom = findall(p -> p[2] < -l / 2 + Δx, eachcol(pc.position))
bc_bottom = PD.VelocityBC(t -> -20, set_bottom, 2)
bc_top = PD.VelocityBC(t -> 20, set_top, 2)

td = PD.VelocityVerlet(2000)

jobname = "test"
path = joinpath("results", jobname)
!ispath(path) && mkpath(path) # create the path if it does not exist
es = PD.ExportSettings(path, 100)

sim = PD.PDSingleBodyAnalysis(name = jobname, pc = pc, mat = mat, precracks = [precrack],
	bcs = [bc_bottom, bc_top], td = td, es = es)

body = PD.init_body(sim.mat, sim.pc)
for precrack in sim.precracks
	PD.define_precrack!(body, precrack)
end

PD.update_thread_cache!(body)
PD.calc_damage!(body)

PD.init_time_discretization!(sim.td, body, sim.mat)

PD.export_vtk(body, sim.es.resultfile_prefix, 0, 0.0)

p = Progress(sim.td.n_steps;
	dt = 1,
	desc = "Time integration... ",
	barlen = 30,
	color = :normal,
	enabled = !PD.is_logging(stderr),
)
Δt½ = 0.5 * sim.td.Δt
begin
    for t in 1:sim.td.n_steps
        time = t * sim.td.Δt
        PD.update_velhalf!(body, Δt½)
        PD.apply_bcs!(body, sim.bcs, time)
        PD.update_disp_and_position!(body, sim.td.Δt)
        PD.compute_forcedensity!(body, sim.mat)
        PD.update_thread_cache!(body)
        PD.calc_damage!(body)
        PD.compute_equation_of_motion!(body, Δt½, sim.mat)
        if mod(t, sim.es.exportfreq) == 0
            PD.export_vtk(body, sim.es.resultfile_prefix, t, time)
        end
        next!(p)
    end
finish!(p)
end
