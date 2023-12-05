function [F,M]=getCoef(Fsection, com)
  # output Cx,Cy,Cz,Cl,Cm,Cn using F_section
  Ftmp = Mtmp = [0 0 0];
  
  x0 = com.x(com.n + 1,:); # center

  for i=1:2*com.n

    Ftmp += Fsection(i,:);
    Mtmp = Mtmp + cross((com.xbv(i,:)-x0), Fsection(i,:));
  end

  # real force = rho * u * (u*gamma)  * (b/2)**2
  # nondim. by 0.5*rho*u**2*real.S
  # real.S = com.S * (b/2)**2

  F = Ftmp/(0.5*com.S);
  M = Mtmp/com.S;
  M(2) = M(2)/(0.5*com.MAC);
  # --> M=[M(1)/(0.5*com.S*2), M(2)/(0.5*com.S*com.MAC), M(3)/(0.5*com.S*2)];
  
endfunction

