function getu=getU(alpha,beta,p,r, phi,windGrad, d_aile, eta_aile, x)

  # z = z/b (normallized)
  # windGrad = [d_u/d_z, d_beta/d_z]

  # aile : aileron deflection angle
  # eta_aile : eta_coordinate where full-moving aileron begin

  omega=[ p 0 r ];

  for i=1:rows(x),
    zeta = x(i,3) + x(i,2)*sin(phi);
    du_wg =  zeta * windGrad(1); # x(3) = dz
    dv_wg = -zeta * windGrad(2)*cos(phi);
    dw_wg =  zeta * windGrad(2)*sin(phi);

    if(x(i,2) >= eta_aile),
      dw_aile = d_aile; # only x>0 side. rad.
    else
      dw_aile = 0.0;
    endif
    
    getu(i,:)=...
     [-cos(alpha)*cos(beta), -cos(alpha)*sin(beta), -sin(alpha)]...
     + [du_wg dv_wg dw_wg  ]...
     + [  0.0   0.0 dw_aile]...
     + cross(x(i,:),omega);
  end
endfunction
 
