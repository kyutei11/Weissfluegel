function dCl = get_dCl(alpha, com)
  beta = p = r = phi = d_aile = 0.0;
  
  [spanwise, total] = getProperty(alpha,beta,p,r, phi, d_aile, com);
  F = total.F_Gamma_bv +total.F_Gamma_di +total.F_Cd_AF + total.F_tv;
  Cl = -F(3)*cos(alpha) + F(1)*sin(alpha);
  
  dCl = Cl - com.targetCl;
endfunction
