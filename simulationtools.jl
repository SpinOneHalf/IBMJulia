using LinearAlgebra
using Interpolations
import Base.Threads.@threads
using Plots
include("utils.jl")
# Parameters
gamma = 1.4
CFL = 0.5

function update(Qnew, dt, dx, dy, nx, ny, Fh, Gh, fb)
        @threads for j = 2 :ny - 1
            for i = 2: nx - 1
                Qnew[i, j, :] = Qnew[i, j, :] - (dt / dx) * (Fh[i, j, :] - Fh[i - 1, j, :]) - (dt / dy) * (
                    Gh[i, j, :] - Gh[i, j - 1, :]) + fb[i, j, :] * dt
            end
        end
    return Qnew
end

function simulation(r, u, v, E, p, c,
                    dx, dy, ny, nx, X0, tMax,dtheta,k=.001)
    # TODO: Redesign around generic force on fluid
    t = 0
    C=c
    while t < tMax
        dt = CFL * minimum([dx, dy]) / maximum(C.+sqrt.(u.*u+v.*v))
        print("time is:")
        println(t)
        Qnew, F, G = generateQFG(r, u, v, E, p)
        Fh, Gh = generateFhGh(Qnew, F, G, c, u, v)

        X0 = X0 + dt / 2 * forcesBody(Qnew, X0, dx, dy, ny, nx)
        FB = forceFluid(dx, dy, nx, ny, X0, dtheta, t,k,u,v)
        Qnew = update(Qnew, dt, dx, dy, nx, ny, Fh, Gh, FB)
        X0 = X0 + dt / 2 * forcesBody(Qnew, X0, dx, dy, ny, nx)
        Qnew[:, 1, :] = Qnew[:, 2, :]
        Qnew[1, :, :] = Qnew[2, :, :]
        Qnew[:, ny - 1, :] = Qnew[:, ny - 2, :]
        Qnew[nx - 1, :, :] = Qnew[nx - 2, :, :]
        # Primitive variables
        r = Qnew[:, :, 1]
        u = Qnew[:, :, 2] ./ r
        v = Qnew[:, :, 3] ./ r
        E = Qnew[:, :, 4]
  
        p = (gamma - 1) * (E - 0.5 * r .* (u .* u + v.* v))
        
        
        C = sqrt.(gamma* p ./ r)
        t = t + dt  # Ardvance  time

    end
    return r,u,v,E,p
end

function generateQFG(r, u, v, E, p)
    nx, ny = size(r)
    Q = zeros((nx, ny, 4))
    F = zeros((nx, ny, 4))
    G = zeros((nx, ny, 4))
    @threads for j =1: ny
        for i = 1: nx
            # Flux calculation
            Q[i, j, 1] = r[i, j]
            Q[i, j, 2] = r[i, j] * u[i, j]
            Q[i, j, 3] = r[i, j] * v[i, j]
            Q[i, j, 4] = E[i, j]

            F[i, j, 1] = r[i, j] * u[i, j]
            F[i, j, 2] = r[i, j] * u[i, j] ^ 2 + p[i, j]
            F[i, j, 3] = r[i, j] * u[i, j] * v[i, j]
            F[i, j, 4] = u[i, j] * (E[i, j] + p[i, j])

            G[i, j, 1] = r[i, j] * v[i, j]
            G[i, j, 2] = r[i, j] * v[i, j] * u[i, j]
            G[i, j, 3] = r[i, j] * v[i, j] ^ 2 + p[i, j]
            G[i, j, 4] = v[i, j] * (E[i, j] + p[i, j])
        end
    end
    return Q, F, G
end

function generateFhGh(Q, F, G, c, u, v)
    nx, ny = size(c)

    ah1 = zeros((nx, ny))
    ah2 = zeros((nx, ny))

    Fh = zeros((nx - 1, ny, 4))
    Gh = zeros((nx, ny - 1, 4))

    @threads for j = 1: ny - 1
        for i =1 :nx - 1
            ah1[i, j] = max(abs(u[i, j]) + c[i, j], abs(u[i + 1, j]) + c[i + 1, j])
            ah2[i, j] = max(abs(v[i, j]) + c[i, j], abs(v[i, j + 1]) + c[i, j + 1])

            Fh[i, j, :] = 0.5 * (F[i + 1, j, :] + F[i, j, :] - ah1[i, j] * (Q[i + 1, j, :] - Q[i, j, :]))
            Gh[i, j, :] = 0.5 * (G[i, j, :] + G[i, j + 1, :] - ah2[i, j] * (Q[i, j + 1, :] - Q[i, j, :]))
        end
    end
    return Fh, Gh
end

function forceFluid(dx, dy, nx, ny, X0, dtheta, t,k,u,v)
    Fk = _curve_force(X0, dtheta, t,k)
    p =size(X0)[1]
    Fb = zeros((nx, ny, 4))
    @threads for i = 1:nx
        for j = 1:ny
            xp = i * dx
            yp = j * dy
            ftemp = zeros((3))
            for l = 1:p
                xtemp = X0[l, :]
                xpt = (floor(xtemp[1] / dx)% nx )* dx
                ypt = (floor(xtemp[2] / dy)% ny) * dy
                w = SixPointDelta((xpt - xp)/(dx))*SixPointDelta((ypt - yp)/dy) / (dx * dy)
                Fl = w * Fk[l, :] * dtheta
                ftemp[1:2] += Fl
                ftemp[3]+=Fl[1]*u[i,j]+Fl[2]*v[i,j]
            Fb[i, j, 2] = ftemp[1]
            Fb[i, j, 3] = ftemp[2]
            Fb[i,j,4]=ftemp[3]
            end
        end
    end
    return Fb
end


function _curve_force(xs, dthe, t, K=.1)
    k = size(xs)[1]
    tempf = zeros((k, 2))
    @inbounds for i = 2:k - 1
        tempf[i, :] = ((xs[i + 1,:] - 2 * xs[i, :] + xs[i - 1,:])./ (dthe ^ 2)) *K
    end
    tempf[1, :] = ((xs[end, :] - 2 * xs[1, :] + xs[2, :])./ (dthe ^ 2)) * K
    tempf[end, :] = ((xs[end-1, :] - 2 * xs[end, :] + xs[1, :])./ (dthe ^2)) * K
    return tempf
end

function forcesBody(Qnew, X0,
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
    vs = zeros((points_n, 2))
    @inbounds for i = 1:points_n
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


