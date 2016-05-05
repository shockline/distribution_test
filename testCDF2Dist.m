n = 100;
x = gamrnd(2,1,n,1);

%Compute the ECDF of x.

x = sort(x);
pEmp = ((1:n)-0.5)' ./ n;

a0 = 1; b0 = 1;
p0Fit = gamcdf(x,a0,b0);
a1 = 2; b1 = 1;
p1Fit = gamcdf(x,a1,b1);
plot([0 1],[0 1],'k--', pEmp,p0Fit,'b+', pEmp,p1Fit,'r+');
xlabel('Empirical Probabilities');
ylabel('(Provisionally) Fitted Gamma Probabilities');
legend({'1:1 Line','a=1, b=1', 'a=2, b=1'}, 'location','southeast');

wgt = 1 ./ sqrt(pEmp.*(1-pEmp));
gammaObj = @(params) sum(wgt.*(gamcdf(x,exp(params(1)),exp(params(2)))-pEmp).^2);
paramHat = fminsearch(gammaObj,[log(a1),log(b1)]);
paramHat = exp(paramHat)