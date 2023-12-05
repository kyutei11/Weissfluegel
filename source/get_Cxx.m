function [dCxyz, dClmn] = get_Cxx(total0, total1, com, d_param01)

 [Cxyz0, Clmn0] = get_Cx(total0, com);
 [Cxyz1, Clmn1] = get_Cx(total1, com);

  dCxyz = (Cxyz1-Cxyz0)/d_param01;
  dClmn = (Clmn1-Clmn0)/d_param01;
  
endfunction
