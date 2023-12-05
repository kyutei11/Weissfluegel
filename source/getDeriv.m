# Get Xac,Cla,Cdi, and lateral coefficient.

# at first set alpha for steady state flight condition, using fzero
alpha0 = 0.0 # init. value for fzero


deltaa=0.01;
deltab=0.01;
deltap=0.01;
deltar=0.01;
delta_windGrad = 0.01;

#get Cla,Xac,Cdi.
beta=0;
p=0;
r=0;
phi = 0;

[cF0,cM0]=getcoeff(n,   0.0,beta,p,r, phi, com);
[cFa,cMa]=getcoeff(n,alpha0,beta,p,r, phi, com);

beta=deltab;
p=0;
r=0;
[cFb0,cMb0]=getcoeff(n,   0.0,beta,p,r, phi, com);
[cFba,cMba]=getcoeff(n,alpha0,beta,p,r, phi, com);

beta=0;
p=deltap;
r=0;
[cFp0,cMp0]=getcoeff(n,   0.0,beta,p,r, phi, com);
[cFpa,cMpa]=getcoeff(n,alpha0,beta,p,r, phi, com);

beta=0;
p=0;
r=deltar;
[cFr0,cMr0]=getcoeff(n,   0.0,beta,p,r, phi, com);
[cFra,cMra]=getcoeff(n,alpha0,beta,p,r, phi, com);

cl=-(cFa(3)-cF0(3));
cla=cl/alpha0
cma=(cMa(2)-cM0(2))/alpha0
xac=0.25-cma/cla
MAC

h=xac-0.25;

#              beta            p              r
Ldlv0=       [cFb0(2)        cFp0(2)        cFr0(2)          # y
              cMb0(1)        cMp0(1)        cMr0(1)          # l
              cMb0(3)        cMp0(3)        cMr0(3) ] ...    # n
 ./([1 1 1]'*[deltab         deltap         deltar]);

Ldlva=  [cFba(2)-cFb0(2) cFpa(2)-cFp0(2) cFra(2)-cFr0(2)     # y
         cMba(1)-cMb0(1) cMpa(1)-cMp0(1) cMra(1)-cMr0(1)     # l
         cMba(3)-cMb0(3) cMpa(3)-cMp0(3) cMra(3)-cMr0(3)]... # n
./([1 1 1]'*[cl*deltab      cl*deltap       cl*deltar]);

Ldlva(4,:)=[cMba(2)-cMb0(2)+cl*h...
            cMpa(2)-cMp0(2)+cl*h...
            cMra(2)-cMr0(2)+cl*h]... # %m?-!
./[cl*deltab      cl*deltap       cl*deltar];

Ldlv0(3,:)=Ldlv0(3,:)+Ldlv0(1,:)*h*MAC/2;
Ldlv0(:,3)=Ldlv0(:,3)+Ldlv0(:,1)*h*MAC;

Ldlva(3,:)=Ldlva(3,:)+Ldlva(1,:)*h*MAC/2;
Ldlva(:,3)=Ldlva(:,3)+Ldlva(:,1)*h*MAC;

Ldlv0,Ldlva
