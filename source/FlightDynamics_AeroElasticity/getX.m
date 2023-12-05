com.n=16;
dth=90/com.n;

th=[-90:dth:90]*pi/180;
eta=sin(th);  # 1..2n+1

th=[-90+dth/2:dth:90]*pi/180;
eta2=sin(th); # 1..2n

for i=1:2*com.n+1,
  xtmp=get_x_c(eta(i),com);
  c  = xtmp(4);
  
  com.x(i,:)=xtmp(1:3);
  com.xtv(i,:)=[com.x(i,1)-c/3, com.x(i,2), com.x(i,3)];

  ##add. v.0.22.10.01
  ea = xtmp(5);
  com.xea(i,:)=[com.x(i,1)-c*(ea-0.25), com.x(i,2), com.x(i,3)]; ## cleared later
  if i>1
	com.dSpar(i-1,:) = com.xea(i,:) - com.xea(i-1,:);
	com.SparDirection(i-1,:) = com.dSpar(i-1,:)/norm(com.dSpar(i-1,:));
  end
end

com.S=0;
com.MAC=0;
clear com.xea

for i=1:2*com.n,
  xtmp=get_x_c(eta2(i),com);
  
  com.xbv(i,:)=xtmp(1:3);
  com.c(i)  = xtmp(4);
  
  com.xbc(i,:)=[com.xbv(i,1)-com.c(i)*0.5, com.xbv(i,2), com.xbv(i,3)];
  ##com.xbc(i,:)=[com.xbv(i,1)-com.c(i)*0.5*(5.73/6.28), com.xbv(i,2), com.xbv(i,3)];

  com.nvect(i,:) = cross(com.x(i,:)-com.xbc(i,:), com.x(i+1,:)-com.xbc(i,:)) / (com.c(i)*0.5); # = dy * nvect
  com.dy(i) = norm(com.nvect(i,:));
  com.nvect(i,:) = com.nvect(i,:)/com.dy(i);
  
  com.S=com.S+com.c(i)*com.dy(i);
  com.MAC=com.MAC + com.c(i)^2*com.dy(i);

  ## added v.0.21.10.01
  ea = xtmp(5);
  com.xea(i,:)=[com.xbv(i,1)-com.c(i)*(ea-0.25), com.xbv(i,2), com.xbv(i,3)];
  com.GJ(i) = xtmp(6);

end
com.AR=4/com.S;
com.MAC=com.MAC/com.S;
com.real.MAC = com.MAC*com.real.b/2;
