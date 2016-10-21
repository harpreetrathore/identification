// Determination of OE parameters as described in Example 6.28 on page 209.
// 6.17

exec('armac1.sci',-1);
b = [0 0 0 0 0 0 0 0.44 0 0.76];
f = [1 0.17]; 
c = 1; d = 1; 
process_oe = armac1(1,b,c,d,f,0.05); 
u = prbs_a(2555,250);
xi = rand(1,2555,'normal');
y = arsimul(process_oe,[u xi]);
z = [y(1:length(u))' u'];
zd = detrend(z,'constant');

// Compute IR for time-delay estimation
exec('cra.sci',-1);
[ir,r,cl_s] = cra(zd);

// Time-delay = 2 samples
// Estimate ARX model (assume known orders)
nf = 1; nb =3; nk = 7;
// tic();
exec('oe.sci',-1);
[thetaN,covfN,nvar,resid] = oe(zd,nb,nf,nk);
// toc()

exec('oe_2.sci',-1);
[thetaN,covfN,resid] = oe_2(zd,[nb,nf,nk]);

// Residual plot
exec('stem.sci',-1);
[cov1,m1] = xcov(resid,24,"coeff");
xset('window',1); 
subplot(2,1,1)
stem(0:24,cov1(25:49)');xgrid();
xtitle('Correlation function of residuals from output y1','lag');
[cov2,m2] = xcov(resid, zd(:,2),24,"coeff");
subplot(2,1,2)
stem(-24:24,cov2');xgrid();
xtitle('Cross corr. function between input u1 and residuals from output y1','lag');

