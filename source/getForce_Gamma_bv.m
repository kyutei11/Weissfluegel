function Force=getForce_Gamma_bv(ub, Gamma, dGamma, com)
  # get forces due to Gamma(trailing vortex); this includes Cdi

  # caution:
  # set aileron deflection angle to 0.0
  # to prevent computational error

  Force = zeros(size(ub));
  
  for i_dest = 1:2*com.n,
    x1 = com.x(i_dest  , :);
    x2 = com.x(i_dest+1, :);
    
    Force(i_dest, :) = cross(ub(i_dest,:), x2-x1) * Gamma(i_dest);
    
  end

endfunction

