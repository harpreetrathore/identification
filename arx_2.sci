
function [theta_bj,opt_err,resid] =  arx_2(varargin)
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
	na = n(1); nb = n(2); //nk = n(3); //nf = n(4);
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
        e = YDATA - _objfun(UDATA,YDATA,p,na,nb,nk);
    endfunction
    tempSum = na+nb
    p0 = linspace(0.1,0.9,tempSum)';
    disp(p0)
    [var,errl] = lsqrsolve(p0,G,size(UDATA,"*"));
    disp(var)
    err = (norm(errl)^2);
    opt_err = err;
	resid = G(var,[]);
    a = 1-poly([var(nb+1:nb+na)]',"q","coeff");
    b = poly([repmat(0,nk,1);var(1:nb)]',"q","coeff");
//    c = poly([1; var(nb+1:nb+nc)]',"q","coeff");
//    d = poly([1; var(nb+nc+1:nb+nc+nd)]',"q","coeff");
//    f = poly([1; var(nb+nd+nc+1:nd+nc+nf+nb)]',"q","coeff");
    p = struct('B',b,'A',a);
    disp('Discrete time BJ model:  y(t) = [B(q)/F(q)]u(t) + [C(q)/D(q)]e(t)')
    theta_bj = p;
    disp(theta_bj)
endfunction

function yhat = _objfun(UDATA,YDATA,x,na,nb,nk)
    x=x(:)
     q = poly(0,'q')
    tempSum = nb+na
    // making polynomials
    b = poly([repmat(0,nk,1);x(1:nb)]',"q","coeff");
    a = 1 - poly([x(nb+1:nb+na)]',"q","coeff")
//    c = poly([1; x(nb+1:nb+nc)]',"q","coeff");
//    d = poly([1; x(nb+nc+1:nb+nc+nd)]',"q","coeff");
//    f = poly([1; x(nb+nd+nc+1:nd+nc+nf+nb)]',"q","coeff");
//    bd = coeff(b*d); cf = coeff(c*f); fc_d = coeff(f*(c-d));
//    if size(bd,"*") == 1 then
//        bd = repmat(0,nb+nd+1,1)
//    end
    aSize = coeff(a);bSize = coeff(b)
    maxDelay = max([length(aSize) length(bSize)])
    yhat = [YDATA(1:maxDelay)]
    for k=maxDelay+1:size(UDATA,"*")
        tempB = 0
        for ii = 2:size(bSize,'*')
            tempB = tempB + bSize(ii)*UDATA(k-ii+1)
        end
        tempA = 0
        for ii = 1:size(aSize,"*")
            tempA = tempA + aSize(ii)*YDATA(k-ii)
        end
//        bdadd = 0
//        for i = 2:size(bd,"*")
//            bdadd = bdadd + bd(i)*UDATA(k-i+1)
//        end
//        fc_dadd = 0
//        for i = 2:size(fc_d,"*")
//            fc_dadd = fc_dadd + fc_d(i)*YDATA(k-i+1)
//        end
//        cfadd = 0
//        for i = 2:size(cf,"*")
//            cfadd = cfadd + cf(i)*yhat(k-i+1)
//        end
        yhat = [yhat; [ tempA+tempB ]];
    end
endfunction
