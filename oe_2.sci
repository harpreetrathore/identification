
function [theta_bj,opt_err,resid] =  oe_2(varargin)
//
	[lhs , rhs] = argn();	
//	if ( rhs < 2 ) then
//			errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while should be 2"), "bj", rhs);
//			error(errmsg)
//	end
//
	z = varargin(1)
//	if ((~size(z,2)==2) & (~size(z,1)==2)) then
//		errmsg = msprintf(gettext("%s: input and output data matrix should be of size (number of data)*2"), "bj");
//		error(errmsg);
//	end
//
//	if (~isreal(z)) then
//		errmsg = msprintf(gettext("%s: input and output data matrix should be a real matrix"), "bj");
//		error(errmsg);
//	end
//
	n = varargin(2)
//	if (size(n,"*")<4| size(n,"*")>5) then
//		errmsg = msprintf(gettext("%s: The order and delay matrix [nb nc nd nf nk] should be of size [4 5]"), "bj");
//		error(errmsg);
//	end
//
//	if (size(find(n<0),"*") | size(find(((n-floor(n))<%eps)== %f))) then
//		errmsg = msprintf(gettext("%s: values of order and delay matrix [nb nc nd nf nk] should be nonnegative integer number "), "bj");
//		error(errmsg);
//	end
//
	nf= n(1); nb = n(2); //nk = n(3); //nf = n(4);
//	
	if (size(n,"*") == 2) then
		nk = 1
	else
		nk = n(3);
	end

    // storing U(k) , y(k) and n data in UDATA,YDATA and NDATA respectively 
    YDATA = z(:,1);
    UDATA = z(:,2);
    NDATA = size(UDATA,"*");
    function e = G(p,m)
        e = YDATA - _objfun(UDATA,YDATA,p,nf,nb,nk);
    endfunction
    tempSum = na+nb
    p0 = linspace(0.1,0.9,tempSum)';
    disp(p0)
    [var,errl] = lsqrsolve(p0,G,size(UDATA,"*"));
    disp(var)
    err = (norm(errl)^2);
    opt_err = err;
	resid = G(var,[]);
    f = poly([1 var(nb+1:nb+nf)]',"q","coeff");
    b = poly([repmat(0,nk,1);var(1:nb)]',"q","coeff");
    p = struct('B',b,'F',f);
    disp('Discrete time model: y(t) = [B(x)/F(x)]u(t) + e(t)');
    theta_bj = p;
    disp(theta_bj)
endfunction

function yhat = _objfun(UDATA,YDATA,x,nf,nb,nk)
    x=x(:)
     q = poly(0,'q')
    tempSum = nb+nf
    // making polynomials
    b = poly([repmat(0,nk,1);x(1:nb)]',"q","coeff");
    f =  -1*poly([x(nb+1:nb+nf)]',"q","coeff")
    fSize = coeff(f);bSize = coeff(b)
    maxDelay = max([length(fSize) length(bSize)])
    yhat = [YDATA(1:maxDelay)]
    for k=maxDelay+1:size(UDATA,"*")
        tempB = 0
        for ii = 2:size(bSize,'*')
            tempB = tempB + bSize(ii)*UDATA(k-ii+1)
        end
        tempF = 0
        for ii = 1:size(fSize,"*")
            tempF = tempF + fSize(ii)*yhat(k-ii)
        end
        yhat = [yhat; [ tempF+tempB ]];
    end
endfunction
