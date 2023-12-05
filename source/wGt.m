function wGt = wGt(x0, x1)
    # output velocity at x0, which is
    # induced by unit vortex from [-inf, x1(2), x1(3)] to [x1(1), x1(2), x1(3)]

    d=cross(x1-x0,x1-x0+[1 0 0]);
    wGt=( 1+((x1-x0)/norm(x1-x0))*[1 0 0]' ) *d /(4*pi*norm(d)^2);

endfunction
