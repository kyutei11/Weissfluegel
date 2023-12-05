# 1-Dim nonlinear solver
# version 0.08.09.04

function [x,itl] = fzero(f,x0,para,tol,itlmax)

    if nargin<5,
        itlmax=2*length(x0);
    endif
    
	if nargin<4,
		tol = 1e-6;
	endif

	f0 = feval(f,x0,para);
	x1 = x0 + 0.01;
	f1 = feval(f,x1,para);
		
	for itl=1:itlmax,
		dx = x1-x0;
		x0 = x1;
		x1 = x1-dx*f1/(f1-f0);
		f0 = f1;
		f1 = feval(f,x1,para);
		
		if norm(f1) < tol,
			x = x1;
			return
		endif
	end

endfunction

