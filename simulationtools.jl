using LinearAlgebra
using Interpolations
import Base.Threads.@threads
using Plots
include("utils.jl")
# Parameters
global gamma = 1.4
global CFL = 0.5

function update!(Qnew, dt, dx, dy, nx, ny, Fh, Gh, fb)
        Qnew[2:nx-1, 2:ny-1, :] = Qnew[2:nx-1, 2:ny-1, :] - 
        (dt / dx) * (Fh[2:nx-1, 2:ny-1, :] - Fh[1:nx-2, 2:ny-1, :]) - 
        (dt / dy) * (Gh[2:nx-1, 2:ny-1, :] - Gh[2:nx-1, 1:ny-2, :]) + fb[2:nx-1,2:ny-1,:] * dt

end


function simulation(r, u, v, E, p, c,
                    dx, dy, ny, nx, X0, tMax,dtheta,k=.001)
    # TODO: Redesign around generic force on fluid
    t = 0
    Q = Array{Float32}(undef,nx, ny, 4)
    F = Array{Float32}(undef,nx, ny, 4)
    G = Array{Float32}(undef,nx, ny, 4)
    ah1 = Array{Float32}(undef,nx, ny)
    ah2 = Array{Float32}(undef,nx, ny)
    Fh = Array{Float32}(undef,nx - 1, ny, 4)
    Gh = Array{Float32}(undef,nx, ny - 1, 4)
    Fb = zeros((nx, ny, 4))
    points_n, _ = size(X0)
    vs = Array{Float32}(undef,points_n, 2)
    ks = size(xs)[1]
    tempf = Array{Float32}(undef,ks, 2)
    while t < tMax
        dt = CFL * minimum([dx, dy]) / maximum(c.+sqrt.(u.*u+v.*v))
        print("time is:")
        println(t)
        generateQFG!(Q,F,G,r, u, v, E, p)
        generateFhGh!(ah1,ah2,Fh,Gh,Q, F, G, c, u, v)

        X0 = X0 + dt / 2 * forcesBody(vs,Q, X0, dx, dy, ny, nx)
        forceFluid!(Fb,dx, dy, nx, ny, X0, dtheta, t,k,u,v,tempf)
        update!(Q, dt, dx, dy, nx, ny, Fh, Gh, Fb)
        X0 = X0 + dt / 2 * forcesBody(vs,Q, X0, dx, dy, ny, nx)
        Q[:, 1, :] = Q[:, 2, :]
        Q[1, :, :] = Q[2, :, :]
        Q[:, ny - 1, :] = Q[:, ny - 2, :]
        Q[nx - 1, :, :] = Q[nx - 2, :, :]
        # Primitive variables
        r = Q[:, :, 1]
        u = Q[:, :, 2] ./ r
        v = Q[:, :, 3] ./ r
        E = Q[:, :, 4]
  
        p = (gamma - 1) * (E - 0.5 * r .* (u .* u + v.* v))
        
        
        c = sqrt.(gamma* p ./ r)
        t = t + dt  # Ardvance  time

    end
    return r,u,v,E,p
end

function generateQFG!(Q,F,G,r, u, v, E, p)
    # Flux calculation
    Q[1: nx, 1:ny, 1] = r[1: nx, 1:ny]
    Q[1: nx, 1:ny, 2] = r[1: nx, 1:ny].* u[1: nx, 1:ny]
    Q[1: nx, 1:ny, 3] = r[1: nx, 1:ny].* v[1: nx, 1:ny]
    Q[1: nx, 1:ny, 4] = E[1: nx, 1:ny]

    F[1: nx, 1:ny, 1] = r[1: nx, 1:ny].* u[1: nx, 1:ny]
    F[1: nx, 1:ny, 2] = r[1: nx, 1:ny].* (u[1: nx, 1:ny].^ 2) + p[1: nx, 1:ny]
    F[1: nx, 1:ny, 3] = r[1: nx, 1:ny].* u[1: nx,1:ny].* v[1: nx, 1:ny]
    F[1: nx, 1:ny, 4] = u[1: nx, 1:ny].* (E[1: nx, 1:ny] + p[1: nx, 1:ny])

    G[1: nx,1:ny, 1] = r[1: nx, 1:ny].* v[1: nx, 1:ny]
    G[1: nx, 1:ny, 2] = r[1: nx, 1:ny].* v[1: nx, 1:ny].* u[1: nx, 1:ny]
    G[1: nx, 1:ny, 3] = r[1: nx, 1:ny].* v[1: nx, 1:ny].^ 2 + p[1: nx, 1:ny]
    G[1: nx, 1:ny, 4] = v[1: nx, 1:ny].* (E[1: nx, 1:ny] + p[1: nx, 1:ny])
end

function generateFhGh!(ah1,ah2,Fh,Gh,Q, F, G, c, u, v)
    ah1[1:nx-1, 1:ny-1] = max.(abs.(u[1:nx-1, 1: ny - 1]) + c[1:nx-1, 1: ny - 1], abs.(u[2:nx, 1: ny - 1]) + c[2:nx, 1: ny - 1])
    ah2[1:nx-1, 1:ny-1] = max.(abs.(v[1:nx-1, 1: ny - 1]) + c[1:nx-1, 1: ny - 1], abs.(v[1:nx-1, 2: ny]) + c[1:nx-1, 2:ny])
    Fh[1:nx-1, 1:ny-1, :] = 0.5 * (F[2:nx, 1:ny-1, :] + F[1:nx-1, 1: ny-1, :] - ah1[1:nx-1, 1:ny-1].* (Q[2:nx, 1:ny-1, :] - Q[1:nx-1,1: ny-1, :]))
    Gh[1:nx-1, 1:ny-1, :] = 0.5 * (G[1:nx-1, 1: ny - 1, :] + G[1:nx-1, 2: ny, :] - ah2[1:nx-1, 1: ny - 1].* (Q[1:nx-1, 2: ny, :] - Q[1:nx-1, 1: ny - 1, :]))
end

function forceFluid!(Fb,dx, dy, nx, ny, X0, dtheta, t,k,u,v,tempf)
    Fk = _curve_force(X0, dtheta, t,k,tempf)
    fill!(Fb,0.0)
    xp=(1:nx)*dx
    yp=(1:ny)*dy
    for l = 1:size(X0)[1]
            xtemp = X0[l, :]
            xpt = xtemp[1]
            ypt = xtemp[2]
            w = SixPointDelta.((xp.-xpt)./(dx))*(SixPointDelta.((yp.-ypt)./dy)')./(dx * dy)
            Fu = w.*(Fk[l, 1] * dtheta)
            Fv = w.*(Fk[l,2]*dtheta)
            Fb[:,:,2] += Fu
            Fb[:,:,3] +=Fv
            Fb[:,:,4]+=Fu.*u+Fv.*v
    end
    return Fb
end


function _curve_force(xs, dthe, t, K=.1,tempf=tempf,)
    k = size(xs)[1]
    tempf[2:k - 1, :] = ((xs[3:k ,:] - 2 * xs[2:k - 1, :] + xs[1:k - 2,:])./ (dthe ^ 2)) *K
    tempf[1, :] = ((xs[end, :] - 2 * xs[1, :] + xs[2, :])./ (dthe ^ 2)) * K
    tempf[end, :] = ((xs[end-1, :] - 2 * xs[end, :] + xs[1, :])./ (dthe ^2)) * K
    return tempf
end

function forcesBody(vs,Qnew, X0,
               dx, dy,
               nx, ny)
    x = (1: nx) * dx
    y = (1: ny) * dy
    r = Qnew[:, :, 1]
    u = Qnew[:, :, 2]./ r
    v = Qnew[:, :, 3]./ r
    ut=linear_interpolation((x,y), u)
    vt=linear_interpolation((x,y),v)
    up = ut.(X0[:,1],X0[:,2])
    vp = vt.(X0[:,1],X0[:,2])
    vs[:, 2] = up
    vs[:, 1] = vp
    return vs
end


