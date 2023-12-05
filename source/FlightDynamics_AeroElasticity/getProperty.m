function [spanwise, total] = getProperty(alpha,beta,p,r, phi, d_aile, alpha0, com)
  # get force at each wing section

  eta_aile = com.eta_aileron;

       com.u =getU(alpha,beta,p,r, phi,com.windGrad, d_aile, eta_aile, com.xbc); # used by fzeros
  spanwise.ub=getU(alpha,beta,p,r, phi,com.windGrad, d_aile, eta_aile, com.xbv);
  spanwise.ut=getU(alpha,beta,p,r, phi,com.windGrad, d_aile, eta_aile, com.xtv);

  ##alpha0 = alpha*ones(2*com.n,1); # init.value for fzeros
  ## alfa0 is given from prev. timestep
  
  #gamma=-com.wnG\un;
  [spanwise.alpha itl] = fzeros('get_dw', alpha0, com, 1.e-8, 2*com.n, false);
  ##plot(spanwise.alpha),grid
  ##pause

  spanwise.Cl = spanwise.Cd = spanwise.Cm = zeros(2*com.n,1);
  spanwise.M_Cm_AF = spanwise.F_Gamma_tv = zeros(2*com.n,3);  
  
  # get gamma from alpha
  for i=1:2*com.n
    eta = com.xbv(i,2);
    [spanwise.Cl(i) spanwise.Cd(i) spanwise.Cm(i)] = get_ClCdCm(eta, spanwise.alpha(i), com);
    un = com.u(i,:) * com.nvect(i,:)';
    ux = com.u(i,1);
    u_2D = norm([ux 0 un]);
    spanwise.Gamma(i) = 0.5*u_2D * spanwise.Cl(i) * com.c(i);
  end
  dGamma=[spanwise.Gamma(1) diff(spanwise.Gamma) -spanwise.Gamma(2*com.n)]; # sign : x-positive direction

  ## added v.0.22.10.01
  [spanwise.th_spar spanwise.T_spar] = getSparTwist(spanwise.Cl, spanwise.Cm, com);
  ##plot(spanwise.th_spar*180/pi)

  #sprintf('calculating CDi + CDp. wait a minutes...')
  
  # get forces due to Gamma and Di
  
  # 0.12.2.1: reset aileron deflection angle
  d_aile = 0.0;
  spanwise.ub=getU(alpha,beta,p,r, phi,com.windGrad, d_aile, eta_aile, com.xbv);

  spanwise.F_Gamma_bv = getForce_Gamma_bv(spanwise.ub, spanwise.Gamma, dGamma, com);
  spanwise.F_Gamma_di = getForce_Gamma_di(             spanwise.Gamma, dGamma, com);

  # add forces due to Cdp
  spanwise.F_Cd_AF = getForce_Cd_AF(spanwise.ub, spanwise.Cd, com);
  
  # ------------------------------------------------------------------------------
  # get aerodynamic coefficients from forces above
  [total.F_Gamma_bv, total.M_Gamma_bv] = getCoef(spanwise.F_Gamma_bv, com);
  [total.F_Gamma_di, total.M_Gamma_di] = getCoef(spanwise.F_Gamma_di, com);
  [total.F_Cd_AF, total.M_Cd_AF]       = getCoef(spanwise.F_Cd_AF   , com);
  
  # add momentum due to cm of airfoil©¢
  total.M_AF = dM = M_all = [0 0 0];

  for i=1:2*com.n
    dx = com.x(i+1,:)-com.x(i,:);
    #dx(1) = 0.0; # axis is in y-z plane
    dM = 0.5*norm(spanwise.ub(i,:))**2 * spanwise.Cm(i)* com.c(i)**2;
    spanwise.M_Cm_AF(i,:) = dM*dx;
    total.M_AF = total.M_AF + dM*dx;
  end
  total.M_AF = total.M_AF / com.S;
  total.M_AF(2) = total.M_AF(2) / (0.5*com.MAC);

  # add coef. due to trailing vortices on wing.
  total.F_tv = total.M_tv = [0 0 0];    # use outside vortex of xbv    
  for i=1:2*com.n
    dF=cross(spanwise.ut(i+1,:) ,dGamma(i+1)*[1 0 0]*norm(com.x(i+1,:)-com.xtv(i+1,:))*3*0.75 );
    spanwise.F_Gamma_tv(i,:) = dF;
    
    # xtv=x-[1 0 0]*c/3;
    total.F_tv = total.F_tv + dF;
    total.M_tv = total.M_tv + cross(com.xtv(i,:),dF);
    
  end
  total.F_tv = total.F_tv / (0.5*com.S);
  total.M_tv = total.M_tv / com.S;
  total.M_tv(2) = total.M_tv(2) / (0.5*com.MAC);

endfunction

