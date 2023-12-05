function Force=getForce_Cd_AF(u, Cd_AF, com)
  # get forces due to Cdp of wing section
  # for the sake of simplicity, downwash angle is not considered (it's 2nd order...)

  Force = zeros(size(u));
  
  # add Cdp
  for i=1:2*com.n
    F=0.5*norm(u(i,:))^2 * Cd_AF(i)*com.c(i)*com.dy(i);
    Force(i,:) = F * u(i,:)/norm(u(i,:)); # sign is positive :: u is usually backward.
  end

endfunction

