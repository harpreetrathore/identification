exec('armac1.sci',-1);
a = 1; b = [0 0 0.6 -0.2];
c = [1 -0.3]; d = [1 -0.5]; 
f = [1 -0.7]; 
process_bj = armac1(a,b,c,d,f,0.05); 
u = prbs_a(2555,250);
xi = rand(1,2555,'normal');
y = arsimul(process_bj,[u xi]);
z = [y(1:length(u))' u'];
zd = detrend(z,'constant');


exec bj.sci;
nb = 2; nc = 1; nd = 1; nf = 1; nk = 2;
tempSum = nb+nc+nd+nf
p0 = linspace(1,2,tempSum)';
[theta_bj,opt_err,resid] =  bj(z,[nb,nc,nd,nf,nk]);
