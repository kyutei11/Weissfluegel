function cdp=get_cdp(cl)
  cd0=0.05;
  cl0=0;
  k=0.2;
  cdp=cd0+k.*(cl-cl0).*(cl-cl0);
end
