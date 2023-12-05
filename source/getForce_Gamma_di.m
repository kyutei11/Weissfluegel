function Force=getForce_Gamma_di(Gamma, dGamma, com)
  # get forces due to Di
  
  # comment : 2012.2.1
  # Even if Gamma is given using deflected aileron
  # (e.g. upwind where aileron has positive AOA), 
  # There are no suction effect due to upwind which causes error of Di, 
  # because 'w' is recomputed using Gamma.

  Force = zeros(length(Gamma), 3);
  
  for i_dest = 1:2*com.n,
    x0 = com.xbv(i_dest,:);
    x1 = com.x(i_dest  , :);
    x2 = com.x(i_dest+1, :);
    
    w = [0 0 0];
    for j_source = 1:2*com.n+1,
      x_source = com.x(j_source,:);
      wGtmp = wGt_inf(x0, x_source); # dest, source
      
      if com.GroundEffect,
        xs_mirror = mirror(x_source, com);
        wGtmp = wGtmp - wGt_inf(x0, xs_mirror);
      end
      
      w = w + 0.5 * wGtmp*dGamma(j_source); # 0.5 : actual effect is half of that of inf. vortex line
    end  #  I got w(i) at section.
    
    Force(i_dest, :) = cross(w, x2-x1) * Gamma(i_dest);
    
  end

endfunction

