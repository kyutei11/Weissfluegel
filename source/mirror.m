function x_mirror=mirror(x, com)
  # used if com.GroundEffect = true
  
  x_mirror = [x(1) x(2) 2*com.h - x(3)];

endfunction
