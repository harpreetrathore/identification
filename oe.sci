// OE Model parameter estimation
// Updated(28-9-16)

///////////////////////////////////////
/////////// ARX Model /////////////////


function [thetaN_arx,covt_arx,nvar,res] = arxc(data,na,nb,nk)
az = max(na,nb+nk-1);
zer = zeros(az,1);
zd = data;
// Zeros appended
zd1(:,1) = [zer; zd(:,1)];
zd1(:,2) = [zer; zd(:,2)];
[r,c] = size(zd1);
t = az+1:r; 
yt = zd1(:,1); ut = zd1(:,2);
yt1 = yt'; ut1 = ut'; // row vector
len1 = length(yt1); 
yt2 = zeros(1,len1-az); ut2 = zeros(1,len1-az); 

// arx(Data,[na nb nk]) 
  for i=1:na
    yt2 = [yt2; -yt1(t-i)];
  end;
  for i=nk:nb+nk-1
    ut2 = [ut2; ut1(t-i)]; 
  end;
[r1,c1] = size(yt2); [r2,c2] = size(ut2);
phit = [yt2(2:r1,:); ut2(2:r2,:)]; 
m1 = phit*phit';
[qm,rm] = qr(m1);
m2 = phit*zd(:,1);
thetaN_arx = inv(rm)*qm'*m2;
// thetaN_arx = inv(m1)*m2;
// thetaN_arx = m1\m2;

[r11,c11] = size(thetaN_arx);
a = thetaN_arx(1:na); b = thetaN_arx(na+1:r11);

// Sum of squared residuals

yhat = phit'*thetaN_arx;
res = zd(:,1) - yhat;
N = length(res);
q = rank(phit);
ssr = res'*res;
sig2 = ssr/(N-q);
nvar = sqrt(sig2);
cov_arx = inv(m1);
covt_arx = diag(cov_arx);
endfunction;
///////////////////////////////
//////////////////////////////

///////////////////////////////////////
////////// Model Display /////////////

function disp_mod(N1,covN1)
len = length(covN1);
B1 = pol2str(N1); 
ind = strindex(B1,['+','-']);  
ind = ind - 1;
if ind~=-1 then B2 = strsplit(B1,ind); 
else B2 = B1;  end;
covB = string(covN1);
  
  if ascii(B2(1)) == 32
  B2 = B2(2:len+1); 
  end; 
  
  B3(1) = ' ';
  for i=1:len
    B3(i) = strsubst(B2(i),'*x','(+-' + covB(i) + ')*x');
  end;

  B4 = B3(1); //disp(15); pause
  
  for i=2:len
  B4 = B4 + ' ' + B3(i);
  end;

disp(B4);
endfunction;
///////////////////////////////////////
///////////////////////////////////////

function [thetaN_oe,covN_oe,nvar,resid] = oe(zd,nb,nf,nk)
exec('deconvol.sci',-1);
[thetaN,covfN,nvar,res] = arxc(zd,nf,nb,nk); 
[r1,c1] = size(thetaN);
yt = zd(:,1); ut = zd(:,2);
m=50;
  if nf==0
    thetaN_oe = thetaN;
    covN_oe = covfN;
    resid = res;
  else
    for k=1:m
      a = thetaN(1:nf);
      b = thetaN(nf+1:r1);  
      A = [1 a']; // Filter
      y = yt(1:length(ut))';
      yf = deconvol(y,A);
      uf = deconvol(ut',A); 
      zf = [yf(1:length(uf))' uf']; 
      zdf = detrend(zf,'constant');
      [thetaNf,covf_a,nvar,resid] = arxc(zdf,nf,nb,nk); 
      thetaN = thetaNf;
      a1 = thetaN(1:nf);
      ba = (norm(a-a1))/norm(a1); 
        if ba<0.005
          break;
        end;
    end;
    thetaN_oe = thetaN;
    covN_oe = covf_a;
  end;
[rt,ct] = size(thetaN_oe);
f_oe = [1 thetaN_oe(1:nf)']; b1 = zeros(1,nk);
b_oe = [b1 thetaN_oe(nf+1:rt)'];
cov_f = covN_oe(1:nf); cov_b = covN_oe(nf+1:rt);

x = poly(0,'x');
  if nf ==0
    disp('Discrete time model: y(t) = B(x)u(t) + e(t)');
  else
    disp('Discrete time model: y(t) = [B(x)/F(x)]u(t) + e(t)');
  end;
  
F = poly(f_oe,'x','coeff');
cov_f1 = [0 cov_f'];
B = poly( b_oe,'x','coeff');
cov_b1 = cov_b';

cov_b1 = abs(cov_b1)
cov_f1 = abs(cov_f1)

b = poly([repmat(0,nk,1);thetaN_oe(1+nf:rt)]',"q","coeff");
f = poly([1; thetaN_oe(1:nf)]',"q","coeff");
bCov = poly([repmat(0,nk,1);cov_b1']',"q","coeff");
fCov = poly([cov_f1]',"q","coeff");
p = struct('B',b,'F',f);
pCov = struct('Bcov',bCov,'Fcov',fCov);
thetaN_oe = struct("value",p,"covariance",pCov);
//disp(thetaN_oe1)

//  if nb==0
//    error('Polynomial B is zero');
//  else
//    disp('B(x) = ');
//    pause
//    disp_mod(B,cov_b);
//  end;
//  
//  if nf~=0
//    disp('F(x) = ');
//    disp_mod(F,cov_f1);
//  end;
    bstr = []
    for ii = nk+1:length(b_oe)
         bstr = bstr + string(b_oe(ii))+" (+-"+string(cov_b1(ii-nk)) +")*q^" + string(ii-1)
         bstr = bstr + "  "
    end
    disp('B(q) = ');
    disp(bstr)
    if f_oe(2) < 0 then
        fstr = "1 "
    else
        fstr = "1 +"
    end
    for ii = 2:length(f_oe)
         fstr = fstr + string(f_oe(ii))+" (+-"+string(cov_f1(ii)) +")*q^" + string(ii-1)
          fstr = fstr + "  "
    end
    disp('F(q) = ');
    disp(fstr)
endfunction;

