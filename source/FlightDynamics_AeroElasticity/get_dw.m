## output difference of washdown between 2-D nonlinear wing section and 3_D weissinger method.
## solution is dw == 0.
## modi. v0.22.10.01 : added spar twist angle theta

## alpha_section -> F,M_spar -> theta -> w at ctl.point

function dw = get_dw(alpha_section, com)
  
  dw = un = w_2D_BC = Gamma_section = ub = zeros(2*com.n,1);

  CL_2D = CM_2D = zeros(2*com.n);

  for i=1:2*com.n
    eta = com.xbv(i,2);
    [CL_2D(i), CD_2D, CM_2D(i)] = get_ClCdCm(eta, alpha_section(i), com); # modi: get also CM
  end

  ##------------------------------------------
  ## modi. v.22.09.29

  th_spar = getSparTwist(CL_2D, CM_2D, com); ## only 1st item th_spar
  
  for i=1:2*com.n

    ## un(i) = com.u(i,:) * com.nvect(i,:)' ;
    ## modi v.0.22.09.29
    ## rotate nvect around elastic axis by theta
    nv = com.nvect(i,:)';

    nv = nv + th_spar(i)*cross(com.SparDirection(i,:)', nv);
    
    un(i) = com.u(i,:) * nv/norm(nv) ;
    ## modi. end

    ux    = com.u(i,1);
    u_2D  = norm([ux 0 un(i)]);  # is considered as a real velocity of 2-D section
    ## in this loop, we simply need washdown. so setting uy=0
    ## beta effect related to dihedral is counted on un(i).
    ## beta effect related to  sweep angle is considered in getCoeff later...

    ## 1. alpha(eta) is given as unknown var.
    ## 2. using alpha(eta), Cl(alpha) is given from airfoil data
    ## 3. using Cl, Gamma_bv is given
    ## 4. Gamma_bv -> w_2D_BC at ctrl point is given by 2D manner
    ##    at this step wing torsion angle is added
    
    ## note : nvect is z(downward)-direction !
	## if CL_2D = 2pi*alpha, B.C. w_2D_BC = 0.0
	## otherwise, B.C. allow non-zero w due to nonlinearity(stall etc.) CL(alpha) property
    w_2D_BC(i) = (CL_2D(i)*com.c(i)/(4*pi*norm(com.xbv(i,:) - com.xbc(i,:)) ) - alpha_section(i)) * u_2D;

    Gamma_section(i) = 0.5 * u_2D * CL_2D(i) * com.c(i);
  end

  ## 5. residual is d_w between (w due to all Gamma) and w_2D_BC at eta
  dw = (com.wnG*Gamma_section + un) - w_2D_BC;

endfunction
