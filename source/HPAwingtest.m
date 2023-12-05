com.real.BL  = [ 0.0 2.25 6.5 10.75 15.0 16.0  ]; # root to tip
com.real.xbv = [ 0.0 0.0  0.0  0.0   0.0 -0.075];
com.real.c   = [ 1.2 1.2  1.1  0.9   0.5  0.3  ];

nBL = length(com.real.BL);

com.real.b  = com.real.BL(nBL) * 2;
com.real.lt = 6.5; # tail moment arm

com.eta_BL = com.real.BL / com.real.BL(nBL); # normalize to eta
com.c_BL   = com.real.c  / com.real.BL(nBL);

com.lt = com.real.lt/(com.real.b);
com.Sv = 0.012;

# load airfoil polar in advance
load '../../XfoilAnalysis/dae11polar0.dat' # these data are computed using XFoil
load '../../XfoilAnalysis/dae11polar1.dat'
load '../../XfoilAnalysis/dae21polar.dat'
load '../../XfoilAnalysis/dae31polar.dat'
load '../../XfoilAnalysis/dae51polar.dat'

com.alpha_polar = dae11polar0(:,1);
com.clmap = [dae11polar0(:,2), dae11polar0(:,2) dae11polar1(:,2) dae21polar(:,2) dae31polar(:,2), dae51polar(:,2)];
com.cdmap = [dae11polar0(:,3), dae11polar0(:,3) dae11polar1(:,3) dae21polar(:,3) dae31polar(:,3), dae51polar(:,3)];
com.cmmap = [dae11polar0(:,5), dae11polar0(:,5) dae11polar1(:,5) dae21polar(:,5) dae31polar(:,5), dae51polar(:,5)];


function xyzc=get_x_c(eta, com)
  
  abseta=abs(eta);
  x=0.0;
  y=eta;
  z=-com.zTip_bhalf*(eta**2.0); # parabolic dihedral
  
  c=interp1(com.eta_BL, com.c_BL, abseta);

  xyzc=[x y z c];

endfunction

function [cl,cd,cm]=get_ClCdCm(eta,alpha_rad, com)

    alpha_deg = alpha_rad * 180.0/pi;
    abseta = abs(eta);
    # 2nd-order linear interp in [alpha, eta] plane    
    cl = interp2(com.eta_BL, com.alpha_polar, com.clmap,  abseta,alpha_deg);
    cd = interp2(com.eta_BL, com.alpha_polar, com.cdmap,  abseta,alpha_deg);
    cm = interp2(com.eta_BL, com.alpha_polar, com.cmmap,  abseta,alpha_deg);
      
endfunction


