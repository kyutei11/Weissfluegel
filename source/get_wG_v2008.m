# Get wGamma_bound.
for j=1:com.n, # source
  x1 = com.x(j+1,:);
  x2 = com.x(j  ,:);
  #j
  for i=1:2*com.n, # destination
    x0 = com.xbc(i,:);

    wGbtmp = wGb(x0,x1,x2);

    wGb_x(i,j) = wGbtmp(1);
    wGb_y(i,j) = wGbtmp(2);
    wGb_z(i,j) = wGbtmp(3);

    # opposite side
    wGb_x(2*com.n-(i-1), 2*com.n-(j-1)) =  wGbtmp(1);
    wGb_y(2*com.n-(i-1), 2*com.n-(j-1)) = -wGbtmp(2);
    wGb_z(2*com.n-(i-1), 2*com.n-(j-1)) =  wGbtmp(3);

    if com.GroundEffect,
        x0 = mirror(x0,com);
        wGbtmp=wGb(x0,x1,x2); # 0:dest, 1..2:source

        wGb_x(i,j) -= wGbtmp(1); # '-' : mirror!
        wGb_y(i,j) -= wGbtmp(2);
        wGb_z(i,j) -= wGbtmp(3);

        # opposite side
        wGb_x(2*com.n-(i-1), 2*com.n-(j-1)) -=  wGbtmp(1);
        wGb_y(2*com.n-(i-1), 2*com.n-(j-1)) -= -wGbtmp(2);
        wGb_z(2*com.n-(i-1), 2*com.n-(j-1)) -=  wGbtmp(3);
    end
  end
end

# Get wGamma_trail.
for j=1:com.n+1, # source
  x1 = com.x(j,:);
  #j
  for i=1:2*com.n, # destination
    x0 = com.xbc(i,:);

    wGttmp = wGt(x0, x1); # 0:dest, 1:source

    wGt_x(i,j) = wGttmp(1);
    wGt_y(i,j) = wGttmp(2);
    wGt_z(i,j) = wGttmp(3);

    # opposite side : note... all of trailing vortex vectors are defined by the same direction
    wGt_x(2*com.n-(i-1), 2*(com.n+1)-j) = -wGttmp(1);
    wGt_y(2*com.n-(i-1), 2*(com.n+1)-j) =  wGttmp(2); # bugfix of sign, on 2008/08
    wGt_z(2*com.n-(i-1), 2*(com.n+1)-j) = -wGttmp(3);

    if com.GroundEffect,
      x0 = mirror(x0, com);
      wGttmp = wGt(x0, x1);
 
      wGt_x(i,j) -= wGttmp(1);
      wGt_y(i,j) -= wGttmp(2);
      wGt_z(i,j) -= wGttmp(3);

      # opposite side
      wGt_x(2*com.n-(i-1), 2*(com.n+1)-j) -= -wGttmp(1);
      wGt_y(2*com.n-(i-1), 2*(com.n+1)-j) -=  wGttmp(2);
      wGt_z(2*com.n-(i-1), 2*(com.n+1)-j) -= -wGttmp(3);
    end

  end
end

for i=1:2*com.n
  for j=1:2*com.n
    wG_x = wGb_x(i,j)+ ( wGt_x(i,j)-wGt_x(i,j+1) );
    wG_y = wGb_y(i,j)+ ( wGt_y(i,j)-wGt_y(i,j+1) );
    wG_z = wGb_z(i,j)+ ( wGt_z(i,j)-wGt_z(i,j+1) );
    
    com.wnG(i,j) = [wG_x wG_y wG_z]*com.nvect(i,:)';
  end
end


