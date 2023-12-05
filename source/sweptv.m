com.real.b = 10.0;

function xyzc=get_x_c(eta, com)
  AR= 5.16;
  taper = 1.0;
  swept = 45.0;
  dihedral = 10.0;
  
  c0 = 4.0/((1+taper)*AR);
  c1 = c0 *taper;
  
  lambda0 = tan(swept*pi/180.0);
  Gamma0  = dihedral *pi/180.0;
  
  abseta=abs(eta);
  x=-abseta*lambda0;
  y=eta*cos(Gamma0);
  z=-abseta*sin(Gamma0);
  c=c0 + abseta*(c0-c1);

  xyzc=[x y z c];

endfunction

function [cl,cd,cm]=get_ClCdCm(eta,alpha, com)
  cl = 2*pi*(alpha+5.0*pi/180);

  cd0=0.05;
  cl0=0;
  k=0.2;
  cd=cd0+k.*(cl-cl0).*(cl-cl0);

  cm = 0.1;
  
endfunction


