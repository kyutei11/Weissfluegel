function [spanwise, com] = getLoading(spanwise, com)

  # get structural load Q,M,T
  spanwise.loading.Q = spanwise.loading.M = spanwise.loading.T = zeros(2*com.n,1); # shear, moment, torque
  
  # left wing (not so quite accurate, because center of wing is fixed...)
  for j=1:com.n # source : tip to root
    dF = spanwise.F_Gamma_bv(j,:)+spanwise.F_Gamma_di(j,:)+spanwise.F_Gamma_tv(j,:)+spanwise.F_Cd_AF(j,:);
    dM = spanwise.M_Cm_AF(j,:);
    x_source = com.xbv(j,:);    
    
    for i=j:com.n # destination : [Q M T] at tip(=i) == zero.
      x_dest = com.xbv(i,:);

      # add moment due to F_Gamma
      dMi = dM + vprd(x_source(:)-x_dest(:), dF);
      dx = com.x(i+1,:)-com.x(i,:); # in later version, x_spar will be defined instead of x
      i_spar = dx/norm(dx);
      j_spar = vprd(dx, com.nvect(i,:)); # normal vector to i_spar
      j_spar = j_spar/norm(j_spar);

      spanwise.loading.Q(i) += -com.nvect(i,:)*dF'; # nvect is along z-direction
      spanwise.loading.M(i) += j_spar *dMi';
      spanwise.loading.T(i) += i_spar *dMi';
    end
  end
  
  # right wing
  for j=2*com.n:-1:com.n+1 # source : tip to root
    dF = spanwise.F_Gamma_bv(j,:)+spanwise.F_Gamma_di(j,:)+spanwise.F_Gamma_tv(j,:)+spanwise.F_Cd_AF(j,:);
    dM = spanwise.M_Cm_AF(j,:);
    x_source = com.xbv(j,:);    
    
    for i=j:-1:com.n+1 # destination : [Q M T] at tip(=i) == zero.
      x_dest = com.xbv(i,:);
      
      # add moment due to F_Gamma
      dMi = dM + vprd(x_source(:)-x_dest(:), dF);
      dx  = com.x(i+1,:)-com.x(i,:); # in later version, x_spar will be defined instead of x
      i_spar = dx/norm(dx);
      j_spar = vprd(dx, com.nvect(i,:)); # normal vector to i_spar
      j_spar = j_spar/norm(j_spar);
      
      spanwise.loading.Q(i) += com.nvect(i,:)*dF'; # positive sign
      spanwise.loading.M(i) += -j_spar *dMi'; # negative sign
      spanwise.loading.T(i) +=  i_spar *dMi';
    end
  end
  
  # Gamma = u*gamma
  # real force = rho*u * Gamma = rho*u**2 * gamma
  
  spanwise.loading.Q = spanwise.loading.Q/(0.5*com.S);
  spanwise.loading.M = spanwise.loading.M/com.S;
  spanwise.loading.T = spanwise.loading.T/(0.5*com.MAC * com.S);

  q = 0.5*com.rho*com.real.airspeed**2;
  com.real.loading.Q = spanwise.loading.Q * q*com.real.S;
  com.real.loading.M = spanwise.loading.M * q*com.real.S*com.real.b;
  com.real.loading.T = spanwise.loading.T * q*com.real.S*com.real.MAC;

endfunction
