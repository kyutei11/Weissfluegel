function [th_spar T_spar]=getSparTwist(CL_2D,CM_2D, com)
  ## output Torque and theta for each spar location

  th_spar = T_spar  = zeros(2*com.n,1);

  halfSpan = com.real.b/2; ## normalised by this

  ## for each L/R side
  for LR=0:1
    for i=1:com.n
      ## outer --> inner
    
      ## get M from outer : positive->pitchup
      ## ---------------------------------------------------------
      ## left/Right
      if LR==0 ## Left
        iLR = i;           ## 1 .. n
        jBegin = iLR+1; jEnd = com.n; ## inner from iLR
      else
        iLR = 2*com.n+1-i; ## 2n..n+1
        jBegin = com.n+1; jEnd = iLR-1;
      end

      qi = 0.5*com.rho*(com.real.airspeed*norm(com.u(iLR,:)))**2.0;

	  cSection  = com.c( iLR)*halfSpan;
	  dySection = com.dy(iLR)*halfSpan;

      ## moment : q*cm*c*cosLam * c*(dy/cosLam) = q*cm*c**2*dy
	  dTsection = qi*CM_2D(iLR)*cSection**2 * dySection;
      T_spar(iLR) = T_spar(iLR) + dTsection;

      ## inner spar includes F too
      ##          v---- nvect is z(downward) direction
      dFsection = - qi*CL_2D(iLR)* cSection * dySection * com.nvect(iLR,:);
      
      for j=jBegin:jEnd
        T_spar(j) = T_spar(j) + dTsection*(com.SparDirection(iLR,:)*com.SparDirection(j,:)');  ## inner prod.

        dx = (com.xea(j,:) - com.xbv(iLR,:))*halfSpan;  ## 0.25c i --> elastic axis j
        dMsection = cross(dFsection, dx);
      
        T_spar(j) = T_spar(j) + com.SparDirection(j,:) * dMsection'; ## inner prod.
      end
    end

    ## get theta : CW in y axis
    ## d_th = (T/GJ)*ds

    ## inner --> outer
    for i=1:com.n
      if LR==0 ## Left
        iLR = com.n-i+1; ## n .. 1
        jBegin = 1; jEnd = iLR-1; # outer from iLR
      else
        iLR = i+com.n; ## n+1 .. 2n
        jBegin = iLR+1; jEnd = 2*com.n;
      end

      ds = norm(com.dSpar(iLR,:))*halfSpan;
      d_th = ds * T_spar(iLR)/com.GJ(iLR);

      th_spar(iLR) = th_spar(iLR) + d_th * 0.5; ## itself : x 0.5
      for j=jBegin:jEnd
        th_spar(j) = th_spar(j) + d_th * (com.SparDirection(iLR,:)*com.SparDirection(j,:)'); ## inner prod.
      end
    end

  end ## L/R

  # real force = rho * u * (u*gamma)  * (b/2)**2
  # nondim. by 0.5*rho*u**2*real.S
  # real.S = com.S * (b/2)**2
  
endfunction

