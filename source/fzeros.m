# multidimentional nonlinear equation solver.
# Koichi Takasaki 04/06/2003

function [x,itl] = fzeros(f,x0,para,tol,itlmax)

	if nargin<5,
		itlmax = 2*length(x0);
	end

	if nargin<4,
		tol = 1e-6;
	end

	xold = x0;
	fold = feval(f,xold,para);
	n = length(x0);
	eps = 1e-4;
	H=zeros(n);
	
	for j=1:n,
		xtmp = x0;
		xtmp(j) += eps;
		fdot = feval(f,xtmp,para) - fold;
		H(:,j) = fdot;
	end
	
	H = eps*inv(H);
	
	x = inf;
	
	for i=1:itlmax,
		dx = real(-H*fold);
		xnew = xold + dx;
		fnew = feval(f,xnew,para);
		if norm(fnew) < tol,
			x = xnew;
			itl=i;
			return
		end
		Hy = H*(fnew-fold);
		H = H + ( ((dx-Hy)*dx') * H)/(dx'*Hy);
		xold = xnew;
		fold = fnew;
	end

endfunction
