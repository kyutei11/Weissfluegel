## modi. v.0.22.09.29 : added elastic axis

# com.real.BL  = [ 0.0  2.25 6.5 10.75 15.0 16.0  ]; # root to tip
# com.real.xbv = [ 0.0  0.0  0.0  0.0   0.0 -0.075];
# com.real.c   = [ 1.2  1.2  1.1  0.9   0.5  0.3  ];

# com.ea       = [ 0.35 0.35 0.35 0.35 0.35 0.35 ]; # elastic axis, %MAC
# com.real.GJ  = [ 15.0 15.0 15.0 12.0   8.0 5.0 ]*1.e3 # Nm2/rad


com.real.BL  = [ 0.0  9.5 13.50 16.0 ]; # root to tip
##com.real.BL  = [ 0.0  9.5  13.25 15.0 ]; # root to tip
##com.real.xbv = [ 0.0  0.0   0.0  0.0  ]; # not used anymore
com.real.c   = [ 1.1  1.1   0.8  0.35 ];

com.washout  = [ 0.0  0.0   0.0   5.0 ] * 1.0*pi/180.0;

com.ea       = [ 0.30 0.30 0.30 0.30 ]; # elastic axis  [x MAC]
com.real.GJ  = [ 15.0 10.0   5.0 0.3 ]*1.e3; # Nm2/rad
##com.real.Iyy = []; later...

nBL = length(com.real.BL);

com.real.b  = com.real.BL(nBL) * 2;
com.real.lt = 6.5; # tail moment arm

com.eta_BL = com.real.BL / com.real.BL(nBL); # normalize to eta
com.c_BL   = com.real.c  / com.real.BL(nBL);

com.lt = com.real.lt/(com.real.b);
com.Sv = 0.012;

# load airfoil polar in advance
# load '../../XfoilAnalysis/dae11polar0.dat' # these data are computed using XFoil
# load '../../XfoilAnalysis/dae11polar1.dat'
# load '../../XfoilAnalysis/dae21polar.dat'
# load '../../XfoilAnalysis/dae31polar.dat'
# load '../../XfoilAnalysis/dae51polar.dat'

# com.alpha_polar = dae11polar0(:,1);
# com.clmap = [dae11polar0(:,2), dae11polar0(:,2) dae11polar1(:,2) dae21polar(:,2) dae31polar(:,2), dae51polar(:,2)];
# com.cdmap = [dae11polar0(:,3), dae11polar0(:,3) dae11polar1(:,3) dae21polar(:,3) dae31polar(:,3), dae51polar(:,3)];
# com.cmmap = [dae11polar0(:,5), dae11polar0(:,5) dae11polar1(:,5) dae21polar(:,5) dae31polar(:,5), dae51polar(:,5)];


function xyzc=get_x_c(eta, com)
  
  abseta=abs(eta);
  ##x=0.0;
  y=eta;
  z=-com.zTip_bhalf*(eta**2.0); # parabolic dihedral

  d_c_075 = 0.75*(com.c_BL(:) - com.c_BL(1));

  x  = interp1(com.eta_BL, 2.0*d_c_075, abseta);
  
  c  = interp1(com.eta_BL, com.c_BL   , abseta);
  ea = interp1(com.eta_BL, com.ea     , abseta);
  GJ = interp1(com.eta_BL, com.real.GJ, abseta);

  xyzc=[x y z c ea GJ]; # added ea and GJ

endfunction

function [Cl,Cd,Cm]=get_ClCdCm(eta,alpha_rad, com)

    # alpha_deg = alpha_rad * 180.0/pi;
    # abseta = abs(eta);
    # # 2nd-order linear interp in [alpha, eta] plane    
    # Cl = interp2(com.eta_BL, com.alpha_polar, com.clmap,  abseta,alpha_deg);
    # Cd = interp2(com.eta_BL, com.alpha_polar, com.cdmap,  abseta,alpha_deg);
    # Cm = interp2(com.eta_BL, com.alpha_polar, com.cmmap,  abseta,alpha_deg);

  wash = interp1(com.eta_BL, com.washout, abs(eta));
  
  Cl = 2*pi*(alpha_rad - wash);
  Cd = 0.0;
  Cm = -0.075;
  
endfunction


