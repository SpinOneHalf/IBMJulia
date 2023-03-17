using LinearAlgebra
using Interpolations
import Base.Threads.@threads
using Plots
include("utils.jl")
# Parameters
global gamma = 1.4
global CFL = 0.5

function update!(Qnew, dt, dx, dy, nx, ny, Fh, Gh, fb)
        @threads for j =2:ny-1
        Qnew[2:nx-1, j, :] = Qnew[2:nx-1, j, :] - 
        (dt / dx) * (Fh[2:nx-1, j, :] - Fh[1:nx-2, j, :]) - 
        (dt / dy) * (Gh[2:nx-1, j, :] - Gh[2:nx-1, j-1, :]) + fb[2:nx-1,j,:] * dt
        end
end

function simulation(r, u, v, E, p, c,
                    dx, dy, ny, nx, X0, tMax,dtheta,k=.001)
    # TODO: Redesign around generic force on fluid
    t = 0
    C=c
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
        dt = CFL * minimum([dx, dy]) / maximum(C.+sqrt.(u.*u+v.*v))
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
        
        
        C = sqrt.(gamma* p ./ r)
        t = t + dt  # Ardvance  time

    end
    return r,u,v,E,p
end

function generateQFG!(Q,F,G,r, u, v, E, p)
    nx, ny = size(r)
    @threads for j = 1:ny
    # Flux calculation
    Q[1: nx, j, 1] = r[1: nx, j]
    Q[1: nx, j, 2] = r[1: nx, j].* u[1: nx, j]
    Q[1: nx, j, 3] = r[1: nx, j].* v[1: nx, j]
    Q[1: nx, j, 4] = E[1: nx, j]

    F[1: nx, j, 1] = r[1: nx, j].* u[1: nx, j]
    F[1: nx, j, 2] = r[1: nx, j].* (u[1: nx, j].^ 2) + p[1: nx, j]
    F[1: nx, j, 3] = r[1: nx, j].* u[1: nx,j].* v[1: nx, j]
    F[1: nx, j, 4] = u[1: nx, j].* (E[1: nx, j] + p[1: nx, j])

    G[1: nx,j, 1] = r[1: nx, j].* v[1: nx, j]
    G[1: nx, j, 2] = r[1: nx, j].* v[1: nx, j].* u[1: nx, j]
    G[1: nx, j, 3] = r[1: nx, j].* v[1: nx, j].^ 2 + p[1: nx, j]
    G[1: nx, j, 4] = v[1: nx, j].* (E[1: nx, j] + p[1: nx, j])
    end
    return Q, F, G
end

function generateFhGh!(ah1,ah2,Fh,Gh,Q, F, G, c, u, v)
    nx, ny = size(c)

          @threads for i =1 :nx - 1
            ah1[i, 1:ny-1] = max(abs.(u[i, 1: ny - 1]) + c[i, 1: ny - 1], abs.(u[i + 1, 1: ny - 1]) + c[i + 1, 1: ny - 1])
            ah2[i, 1:ny-1] = max(abs.(v[i, 1: ny - 1]) + c[i, 1: ny - 1], abs.(v[i, 2: ny]) + c[i, 2:ny])

            Fh[i, 1:ny-1, :] = 0.5 * (F[i + 1, 1:ny-1, :] + F[i, 1: ny-1, :] - ah1[i, 1:ny-1].* (Q[i + 1, 1:ny-1, :] - Q[i,1: ny-1, :]))
            Gh[i, 1:ny-1, :] = 0.5 * (G[i, 1: ny - 1, :] + G[i, 2: ny, :] - ah2[i, 1: ny - 1].* (Q[i, 2: ny, :] - Q[i, 1: ny - 1, :]))
        end
    return Fh, Gh
end

function forceFluid!(Fb,dx, dy, nx, ny, X0, dtheta, t,k,u,v,tempf)
    Fk = _curve_force(X0, dtheta, t,k,tempf)
    p =size(X0)[1]
    fill!(Fb,0.0)
    ftemp = zeros((3))
     @threads for i = 1:nx
        ftemp = zeros((3))
        for j = 1:ny
            xp = i * dx
            yp = j * dy
            
            for l = 1:p
                xtemp = X0[l, :]
                xpt = (floor(xtemp[1] / dx)% nx )* dx
                ypt = (floor(xtemp[2] / dy)% ny) * dy
                w = SixPointDelta((xpt - xp)/(dx))*SixPointDelta((ypt - yp)/dy) / (dx * dy)
                Fl = w * Fk[l, :] * dtheta
                ftemp[1:2] += Fl
                ftemp[3]+=Fl[1]*u[i,j]+Fl[2]*v[i,j]
            Fb[i, j, 2:4] .= ftemp
            fill!(ftemp,0.0)
        end
        end
    end
    return Fb
end


function _curve_force(xs, dthe, t, K,tempf)
    k = size(xs)[1]
    @threads for i = 2:k - 1
        tempf[i, :] = ((xs[i + 1,:] - 2 * xs[i, :] + xs[i - 1,:])./ (dthe ^ 2)) *K
    end
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
    points_n, _ = size(X0)
    ut=linear_interpolation((x,y), u)
    vt=linear_interpolation((x,y),v)
    @threads for i = 1:points_n
        temp_x = X0[i, :]
        xp = (floor(temp_x[1] / dx)% nx) * dx
        yp = (floor(temp_x[2] / dy)% ny) * dy
        up = ut(xp,yp)
        vp = vt(xp,yp)
        vs[i, 2] = up
        vs[i, 1] = vp
    end
    return vs
end


