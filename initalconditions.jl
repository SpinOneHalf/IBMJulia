
function generate_circle(r, n,d=.5)
    xs = zeros((n, 2))
    dtheta = (2 * pi) / n
    for i = 1:n
        xs[i, :] = [cos(dtheta * i) * r, sin(dtheta * i) * r].+d
    end
    return xs, dtheta
end
function  static_case(nx=100, ny=100)
    # Initialization
    c = zeros((nx, ny))
    r = zeros((nx, ny))
    u = zeros((nx, ny))
    v = zeros((nx, ny))
    E = zeros((nx, ny))
    p = zeros((nx, ny))
    dx = 1 / (nx - 1)
    dy = 1 / (ny - 1)
    # Creating mesh
    for i = 1:nx
        for j = 1:ny
            p[i, j] = .5
            r[i, j] = .5323
            u[i, j] = 0
            v[i, j] = 0
            E[i, j] = p[i, j] / (gamma - 1) + 0.5 * r[i, j] * (u[i, j] ^2 + v[i, j]^2)
            c[i, j] = sqrt(gamma * p[i, j] / r[i, j])
        end
    end
    return r, u, v, p, E, c, dx, dy
end
