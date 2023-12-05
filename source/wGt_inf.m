function wGt_inf = wGt_inf(x0, x1)
    # output velocity at x0, which is
    # induced by unit vortex from [-inf, x1(2), x1(3)] to [inf, x1(2), x1(3)]

    d=cross(x1-x0, x1-x0 +[1 0 0]);
    wGt_inf= d /(2*pi*norm(d)^2);
endfunction
