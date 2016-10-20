// ARX model parameter estimation
// Computes Covariance matrix
// Computes Noisevariance of a process
// Updated(28-9-16)

function [thetaN,covt,nvar,res] = arx(data,na,nb,nk)
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
[r1_a,c1_a] = size(yt2); [r2_a,c2_a] = size(ut2);
phit = [yt2(2:r1_a,:); ut2(2:r2_a,:)];
m1 = phit*phit';
[qm,rm] = qr(m1);
m2 = phit*zd(:,1);
thetaN = inv(rm)*qm'*m2;
// thetaN = inv(m1)*m2;
// thetaN = m1\m2;

[r1,c1] = size(thetaN);
a = thetaN(1:na); b = thetaN(na+1:r1);

// Sum of squared residuals
yhat = phit'*thetaN;
res = zd(:,1) - yhat;
N = length(res);
q = rank(phit);
ssr = res'*res;
sig2 = ssr/(N-q);
nvar = sqrt(sig2);
cov = inv(m1);
covt = diag(cov);

a = thetaN(1:na); b = [repmat(0,nk,1);thetaN(na+1:r1)]; 
cova = covt(1:na); covb = covt(na+1:r1); 
x = poly(0,'x');
disp('Discrete time model: A(x)y(t) = B(x)u(t) + e(t)');
A = poly([1 a'],'x','coeff');
cov_a1 = [0 cova'];
b1 = zeros(1,nk);
B = poly([b1 b'],'x','coeff');
cov_b1 = covb';

cov_b1 = abs(cov_b1)
cov_a1 = abs(cov_a1)
bpol = poly([repmat(0,nk,1);thetaN(1+na:na+nb)]',"q","coeff");
apol = poly([1; thetaN(1:na)]',"q","coeff");
bCov = poly([repmat(0,nk,1);cov_b1']',"q","coeff");
aCov = poly([cov_a1]',"q","coeff");
p = struct('B',bpol,'F',apol);
pCov = struct('Bcov',bCov,'Acov',aCov);
thetaN = struct("value",p,"covariance",pCov);
    bstr = []
    for ii = nk+1:length(b)
         bstr = bstr + string(b(ii))+" (+-"+string(cov_b1(ii-nk)) +")*q^" + string(ii)
         bstr = bstr + "  "
    end
    disp('B(q) = ');
    disp(bstr)
    if a(1) < 0 then
        astr = "1 "
    else
        astr = "1 +"
    end
    for ii = 1:length(a)
         astr = astr + string(a(ii))+" (+-"+string(cov_a1(ii+1)) +")*q^" + string(ii)
          astr = astr + "  "
    end
    disp('A(q) = ');
    disp(astr)

endfunction;

