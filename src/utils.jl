function put_pts_on_circle(xc, yc, r, n)
    t = range(0, 2π, length=n+1)
    pts = zeros(n, 2)
    for i in 1:n
        pts[i, 1] = xc + r*cos(t[i])
        pts[i, 2] = yc + r*sin(t[i])
    end
    return pts
end

function put_pts_on_ellipse(xc, yc, rx, ry, n)
    t = range(0, 2π, length=n+1)
    pts = zeros(n, 2)
    for i in 1:n
        pts[i, 1] = xc + rx*cos(t[i])
        pts[i, 2] = yc + ry*sin(t[i])
    end
    return pts
end

function reparametrization(path, h)
    dp = path .- circshift(path, (1, 0))
    dp[1, :] .= 0
    dl = sqrt.(sum(dp.^2, dims=2))
    lp = vec(cumsum(dl, dims=1))
    total_len = lp[end]
    lp ./= total_len
    npath = round(Int, total_len / h)
    g1 = range(0, 1, length=npath)
    itp_x = LinearInterpolation(lp, path[:, 1])
    itp_y = LinearInterpolation(lp, path[:, 2])
    path_x = itp_x.(g1)
    path_y = itp_y.(g1)
    return [path_x path_y]
end
