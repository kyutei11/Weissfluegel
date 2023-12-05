# Weissfluegel : 3-D nonlinear stability derivative analysis software
# Version 1.0, Copyright (C) 2004 Koichi Takasaki
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# This program is distributed WITHOUT ANY WARRANTY;
# without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#
# See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#
# E-mail : koichi.takasaki@gmail.com

clear

# load wing geometry and 2-D polar file

HPAwingtest
#HPAwingtest_linGamma

com.GroundEffect = true;
#com.h = 0.1; # -z/0.5b, for ground effect calc.
com.h = 0.3125; # -z/0.5b, for ground effect calc. z = 5.0m.

com.windGrad = [0.000 0.000]; # windgrad_x, windgrad_y

#--------------------------------------------------
com.rho           =  1.2; # kg/m3
com.real.airspeed =  7.0; # m/s
com.real.weight   = 95.0; # kg

#--------------------------------------------------
# dihedral angle (delta_z at tip)
#com.zTip_bhalf = 0.0000;
com.zTip_bhalf = 0.0625;
#com.zTip_bhalf = 0.125;

#com.zTip_bhalf = 0.1468; # for HPAwingtest_linGamma equivalent Clbeta for 0.125(parabolic)
#com.zTip_bhalf = 0.072064; # for HPAwingtest_linGamma equivalent Clbeta for 0.0625(parabolic)

#--------------------------------------------------
# aileron position ( only eta>0 side )
com.eta_aileron = 15000/16000;

#--------------------------------------------------
getX   # get wing geometory
get_wG # get influece matrices

com.real.S  = com.S*(0.5*com.real.b)**2;
com.targetCl = com.real.weight*9.8/(0.5*com.rho * com.real.airspeed**2* com.real.S);

alpha = beta = p = r = phi = 0.0;

[target_alpha, itl] = get_alphaTarget(com.targetCl, 0.0*pi/180, com);

alpha = target_alpha; # for later computation

# step:
# set alpha or beta or p or r etc. += 0.01 etc. ; if += 0.0, save as 00.dat
# call function below
# save com, spanwise and total

#alpha  = target_alpha + 0.0;
#beta   = 0.0;
#p      = 0.0;
#r      = 0.0;
#phi    = 0.0;
#d_aile = 0.0;
#
#[spanwise, total] = getProperty(alpha,beta,p,r, phi, d_aile, com);

###getDeriv # get stability derivative : not used anymore
