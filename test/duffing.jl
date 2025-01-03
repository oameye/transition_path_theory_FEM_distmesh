
using Test
using FEMTPT
using Random
using Plots

# parameters for the Duffing oscillator with noise:
# X'' + gamma*X' + X^3 - X = \sqrt{2\gamma\beta^{-1}}\eta_t

xa = -1.0
ya = 0.0
xb = 1.0
yb = 0.0

rx = 0.3
ry = 0.4
beta = 20.0
Hbdry = 0.5

h0 = 0.1
gamma = 0.5

function Hamiltonian(x, y)
    0.5 .* y.^2 .+ 0.25 .* x.^4 .- 0.5 .* x.^2
end

function KE(x)
    0.5 .* (x[:,2].^2)
end

function divfree(x, y)
    f1 = y
    f2 = .- x.^3 .+ x
    return f1, f2
end

function drift(x, y)
    f1 = y
    f2 = .- gamma .* y .- x.^3 .+ x
    return f1, f2
end

nx, ny = 41, 41
nxy = nx * ny
xmin, xmax = -2.0, 2.0
ymin, ymax = -2.0, 2.0

x1 = range(xmin, xmax, length=nx)
y1 = range(ymin, ymax, length=ny)

x_grid = [xx for yy in y1, xx in x1]
y_grid = [yy for yy in y1, xx in x1]

drift1, drift2 = drift(x_grid, y_grid)
dnorm = sqrt.(drift1.^2 .+ drift2.^2 .+ 1e-12)
Hgrid = Hamiltonian(x_grid, y_grid)

# generate a trajectory

Ntraj = 1_000_000
h = 1e-3 # time step
sqh = sqrt(h * 2 * gamma / beta)

px = zeros(Ntraj)
py = zeros(Ntraj)
w = sqh * randn(Ntraj - 1)
# initial point
px[1] = xa
py[1] = ya
for k in 1:Ntraj-1
    drift1, drift2 = drift(px[k], py[k])
    px[k+1] = px[k] + drift1 * h
    py[k+1] = py[k] + drift2 * h + w[k]
end

Na = round(Int, π * (rx + ry) / h0) # the number of points on the A-circle
Nb = round(Int, π * (rx + ry) / h0) # the number of points on the B-circle

ptsA = put_pts_on_ellipse(xa, ya, rx, ry, Na)
ptsB = put_pts_on_ellipse(xb, yb, rx, ry, Nb)

using Contour
cont = Contour.contour(x1, y1, Hgrid, 0.5)
yc, xc =coordinates(lines(cont)[1])
p_outer = [xc yc]

pts_outer = reparametrization(p_outer,h0);
Nouter = size(pts_outer, 1)
Nfix = Na+Nb+Nouter


# @testset "Duffing" begin
# input data for triangulation
bbox = [xmin, xmax, ymin, ymax]
pfix = zeros(Nfix, 2)
pfix[1:Na, :] .= ptsA
pfix[Na+1:Na+Nb, :] .= ptsB
pfix[Na+Nb+1:Nfix, :] .= pts_outer

dA = dellipse(pfix, xa, ya, rx, ry)
function dfunc(p)
    d0 = Hamiltonian(p[:, 1], p[:, 2])
    dA = dellipse(p, xa, ya, rx, ry)
    dB = dellipse(p, xb, yb, rx, ry)
    d = ddiff(d0 .- Hbdry, FEMTPT.dunion(dA, dB))
    return d
end

pts, tri = distmesh2D(dfunc, FEMTPT.huniform, h0, bbox, pfix)

tri[1,2]
fig = plot(aspect_ratio=:equal)
for i in 1:size(tri,1)
    plot!([pts[tri[i,j],1] for j in [1,2,3,1]],
          [pts[tri[i,j],2] for j in [1,2,3,1]],
          color=:black, linewidth=0.1, legend=false)
end
fig


# FEM_committor_solver_irreversible(pts,tri,Aind,Bind,Fpts,drift,divergence,eps05)
NAind,Aind = find_ABbdry_pts_ellipse(pts,xa,ya,rx,ry,h0) # find mesh points on \partial A
NBind,Bind = find_ABbdry_pts_ellipse(pts,xb,yb,rx,ry,h0) # find mesh points on \partial B


q = FEMTPT.FEM_committor_solver_Langevin(pts,tri,Aind,Bind,KE,divfree,beta,gamma)

# TPTdata = np.concatenate((pts,np.reshape(q,(Npts,1))),axis = 1)

qmin = minimum(q)
qmax = maximum(q)

using TriplotRecipes

plot(aspect_ratio=:equal,size=(800,720))
tripcolor!(pts[:,1], pts[:,2], q, tri',color=:magma)
# trimesh!(pts[:,1], pts[:,2], tri',fillalpha=0.0,linecolor=:white, lw=0.1)

function divfree1(x,y)
    f1,f2 = divfree(x,y)
    return -f1,-f2
end

qminus = FEMTPT.FEM_committor_solver_Langevin(pts,tri,Bind,Aind,KE,divfree1,beta,gamma)

# Triangulate sets A and B.
# This is necessary for finding the normalization constant Z for the invariant density.

# input data for triangulation
# bbox = [xmin,xmax,ymin,ymax]

function dfuncA(p)
    return dellipse(p, xa, ya, rx, ry)
end

function dfuncB(p)
    return dellipse(p, xb, yb, rx, ry)
end

bboxA = [xa-rx, xa+rx, ya-ry, ya+ry]
pts_Amesh, tri_Amesh = distmesh2D(dfuncA, FEMTPT.huniform, h0, bboxA, ptsA)
bboxB = [xb-rx, xb+rx, yb-ry, yb+ry]
pts_Bmesh, tri_Bmesh = distmesh2D(dfuncB, FEMTPT.huniform, h0, bboxB, ptsB)

Npts_Amesh = size(pts_Amesh, 1)
Ntri_Amesh = size(tri_Amesh, 1)
println("Npts = ", Npts_Amesh, " Ntri = ", Ntri_Amesh)

Npts_Bmesh = size(pts_Bmesh, 1)
Ntri_Bmesh = size(tri_Bmesh, 1)
println("Npts = ", Npts_Bmesh, " Ntri = ", Ntri_Bmesh)

function fpot(x)
    return Hamiltonian(x[:, 1], x[:, 2])
end

Z = FEMTPT.invariant_pdf(pts, tri, pts_Amesh, tri_Amesh, pts_Bmesh, tri_Bmesh, fpot, beta)

Rcurrent, Rrate =  FEMTPT.reactive_current_transition_rate_Langevin(pts, tri, fpot, divfree, beta, gamma, q, qminus, Z)

ARcurrent = vec(sqrt.(sum(Rcurrent.^2, dims=2)))
ARCmax = maximum(ARcurrent)

tripcolor(pts[:,1], pts[:,2], ARcurrent, tri',color=:viridis)

c = [ARcurrent ARcurrent]'
quiver(pts[:,1], pts[:,2], quiver=(Rcurrent[:,1], Rcurrent[:,2]),line_z=repeat([c...], inner=2), c=:thermal)

prob_reactive = FEMTPT.probability_reactive_Langevin(pts,tri,fpot,beta,q,qminus,Z)

prob_lastA = FEMTPT.probability_last_A_Langevin(pts,tri,pts_Amesh,tri_Amesh,fpot,beta,qminus,Z)
