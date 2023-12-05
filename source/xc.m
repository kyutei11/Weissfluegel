function xc=xc(eta)
  aeta=abs(eta);
  
  etaref=106.05/150;
  c0=60/150;
  c1=45/150;
  c2=30/150;
  lm0=(60-45)/110*0.75;
  lm1=(45-30)/40 *0.75;
  gm0=10*pi/180;
  gm1=30*pi/180;

  if aeta<etaref,
    x=-aeta*lm0;
    y=  eta*cos(gm0);
    z=-aeta*sin(gm0);
    c=c0+(aeta/etaref)*(c1-c0);
  else
    x=-etaref*lm0-(aeta-etaref)*lm1;
    y=sign(eta)*etaref*cos(gm0)...
      +(eta-sign(eta)*etaref)*cos(gm1);
    z=-etaref*sin(gm0)-(aeta-etaref)*sin(gm1);
    c=c1+((aeta-etaref)/(1-etaref))*(c2-c1);
  end
  xc=[x y z c];
endfunction
