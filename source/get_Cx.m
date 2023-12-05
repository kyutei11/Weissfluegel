function [Cxyz, Clmn] = get_Cx(total, com)
  
  Cxyz = total.F_Gamma_bv + total.F_Gamma_di + total.F_Cd_AF + total.F_tv;
  Clmn = total.M_Gamma_bv + total.M_Gamma_di + total.M_Cd_AF + total.M_tv + total.M_AF;
  
endfunction