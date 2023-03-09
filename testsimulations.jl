include("simulationtools.jl")
include("plottingtools.jl")
include("initalconditions.jl")
nx = 150
ny = 150
r0, u0, v0, p0, E0, c0, dx, dy = static_case(nx, ny)
xs, dtheta = generate_circle(.2, 500)

r,u,v,E,p = simulation(r0, u0, v0, E0, p0, c0,
                    dx, dy, ny, nx, xs, .2,dtheta)
heatmap(p)
print("DONE")