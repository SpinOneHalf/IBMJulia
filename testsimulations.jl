include("simulationtools.jl")
include("initalconditions.jl")
using Plots
nx = 250
ny = 250
r0, u0, v0, p0, E0, c0, dx, dy = static_case(nx, ny)
xs, dtheta = generate_circle(.1, 100)

r,u,v,E,p = simulation(r0, u0, v0, E0, p0, c0,
                    dx, dy, ny, nx, xs, .3,dtheta)
heatmap(p)
print("DONE")