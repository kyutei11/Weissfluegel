function com = getCGz_wing(spanwise, com)

  # CG point along is z-axis of wing
  # is approximated by integ(M*z)dy / integ(M)dy
  
  # type '[spanwise, com] = getLoading(spanwise, com)'
  # before using this function
  
  mz_dy = 0;
  m_dy = 0;

  for i=1:com.n
    mz_dy = mz_dy + spanwise.loading.M(i)*com.xbc(i,3) * com.dy(i);
    m_dy  =  m_dy + spanwise.loading.M(i)              * com.dy(i);
  end  

  com.CGz_wing = mz_dy/m_dy;

endfunction
