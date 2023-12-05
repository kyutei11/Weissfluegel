
com.xzc_section = zeros(4, 5) # [eta,x,z,c] , sectionNo

com.b = 30.0
com.eta_BL(:) = [0, 0.5, com.b] / com.b
com.xzc_section(2,:) = []
# z is not defined, because of parabolic dihedral


function xzc=get_x_c(eta, com)
  span = 30.0
  
  aeta=abs(eta);
  x=0;
  y=eta;
  z=0;
  c=0.2*sqrt(1-aeta^2);

  xyzc=[x y z c];

endfunction

function [cl,cd,cm]=get_clcdcm(eta,alpha)
  cl = 2*pi*alpha

  cd0=0.05;
  cl0=0;
  k=0.2;
  cdp=cd0+k.*(cl-cl0).*(cl-cl0);

  cd = 0.0

  cm = 0.0
endfunction