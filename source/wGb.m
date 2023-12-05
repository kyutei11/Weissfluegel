function wGb= wGb(x0, x1, x2)
    # output velocity at x0, which is
    # induced by unit vortex from x1 to x2

    d=cross(x1-x0,x2-x0)/norm(x2-x1);
    wGb=( ((x1-x0)/norm(x1-x0)-(x2-x0)/norm(x2-x0))...
             *((x2-x1)/norm(x2-x1))' )...
             *d /(4*pi*norm(d)^2);

endfunction
