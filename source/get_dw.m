# output difference of washdown between 2-D nonlinear wing section and 3_D weissinger method.
# solution is dw == 0.

function dw = get_dw(alpha_section, com)
  dw = un = w_2D = Gamma_section =   zeros(2*com.n,1);
  
  for i=1:2*com.n
    eta = com.xbv(i,2);
    CL_2D = get_ClCdCm(eta, alpha_section(i), com); # 1st item
    [eta alpha_section(i) CL_2D];
    
    un(i) = com.u(i,:) * com.nvect(i,:)' ;
    ux    = com.u(i,1);
    u_2D  = norm([ux 0 un(i)]);  # is considered as a real velocity of 2-D section
    
    # note : nvect is z-direction!
    w_2D(i) = (CL_2D*com.c(i)/(4*pi*norm(com.xbv(i,:) - com.xbc(i,:)) ) - alpha_section(i)) * u_2D;
    
    Gamma_section(i) = 0.5 * u_2D * CL_2D * com.c(i);
  end
  
  dw = (com.wnG*Gamma_section + un) - w_2D;

endfunction
