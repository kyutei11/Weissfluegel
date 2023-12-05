## Weissfluegel : 3-D nonlinear stability derivative analysis software
## Version 1.0, Copyright (C) 2004 Koichi Takasaki
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 2
## of the License, or (at your option) any later version.
## This program is distributed WITHOUT ANY WARRANTY;
## without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
##
## See the GNU General Public License for more details.
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
##
## E-mail : koichi.takasaki@gmail.com
##
## ver. 22.09.28 : modified for FlightDynamics and Aeroelasticity coupled analysis

clear

## load wing geometry and 2-D polar file

HPAwingtest
##HPAwingtest_linGamma

com.windGrad = [0.0 0.0];

com.HalfModel    = false;

com.GroundEffect = false;
com.h = 0.3125; ## -z/0.5b, for ground effect calc. z = 5.0m.

com.eta_aileron = 1.0; ## no aileron

##com.zTip_bhalf = 0.125;
com.zTip_bhalf = 0.;

##--------------------------------------------------
com.rho           =  1.2; ## kg/m3
com.real.airspeed = 11.11; ## m/s
com.real.weight   = 95.0; ## kg

##--------------------------------------------------
getX   ## get wing geometory
get_wG ## get influece matrices

com.real.S  = com.S*(0.5*com.real.b)**2;

## init. flight traj.

alpha = 8.0*pi/180; ## about 8.5m/s, little dive
d_aile = 0.0;
beta = p = r = phi = 0.0;

alpha0 = alpha*ones(2*com.n,1); # init.value for fzeros

[spanwise, total] = getProperty(alpha,beta,p,r, phi, d_aile, alpha0, com);

## step:
## set alpha or beta or p or r etc. += 0.01 etc. ; if += 0.0, save as 00.dat
## call function below
## save com, spanwise and total

##alpha  = target_alpha + 0.0;
##beta   = 0.0;
##p      = 0.0;
##r      = 0.0;
##phi    = 0.0;
##d_aile = 0.0;
##

## for t=0:0.1:10.0
								#   [spanwise, total] = getProperty(alpha,beta,p,r, phi, d_aile, alpha0, com);
  

## end
