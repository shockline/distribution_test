function [pEmp, names, pFit]= CDFFit(data,method,maxPointPloted)
if nargin == 1
    method = 'L2';
end

if nargin == 2
    maxPointPloted = 80;
end

switch method
    case 'L1Sum'
        measure = @(params) sum(abs(params));
    case 'L1Max'
        measure = @(params) max(abs(params));
    case 'L2Sum'
        measure = @(params) sum((params).^2);
    otherwise
        measure = @(params) sum((params).^2);
end

x = data(:,1);
pEmp = data(:,2);
names = {'Gamma', 'InverseGaussian', 'Lognormal', 'Rician', ...
    'Burr', 'BirnbaumSaunders', 'Rayleigh', 'Weibull', ...
    'Noncentral Chi-square', 'Noncentral F', 'F', 'Nakagami', ...
    'Logistic', 'LogLogistic',  'Levy'};

% names = {'Gamma', 'InverseGaussian', 'ExtremeValue', 'Lognormal', 'Rician', ...
%     'Burr', 'BirnbaumSaunders', 'Rayleigh', 'Weibull', ...
%     'F', 'Nakagami', ...
%     'Logistic', 'LogLogistic',  'Levy'};

% names = {'Beta'};

% names = {'ExtremeValue'};


legendList = {'b+', 'r+', 'black+', 'm+', 'bpentagram', 'rpentagram', ...
    'blackpentagram', 'mpentagram', 'b*', 'r*', 'black*', 'm*', 'bo', 'ro','blacko'};
pFit = zeros(size(x,1),size(names,2));
plot([0 1],[0 1],'k--')
hold on;
for i=1:size(names,2)
    fprintf('Now processing: %s\n', names{i});
    switch  names{i}        
        case 'Beta'
            [a1, b1] = estimateBetaParameters(x,pEmp,measure);
            pFit(:,i) = cdf('Beta',x,a1,b1)';        
            fprintf('Beta:\ta=%f\tb=%f\n',a1,b1)
        case 'Gamma'
            [k1, theta] = estimateGammaParameters(x,pEmp,measure);
            pFit(:,i) = cdf('Gamma',x,k1,theta)';
        case  'InverseGaussian'
            [miu1, lambda] = estimateIGParameters(x,pEmp,measure);
            pFit(:,i) = cdf('InverseGaussian',x,miu1,lambda)';
            ig_mode = miu1*((1+9*miu1^2/4/lambda^2)^0.5-3*miu1/2/lambda);
            fprintf('IG:\tmiu=%f\tlambda=%f\tmode=%f\n',miu1,lambda,ig_mode)
        case  'ExtremeValue'
            [miu6, sigma6] = estimateExtremeValueParameters(x,pEmp,measure);
            pFit(:,i) = evcdf(x,miu6,sigma6)';
        case 'Lognormal'
            [miu2, sigma1] = estimateLognParameters(x,pEmp,measure);
            pFit(:,i) = cdf('Lognormal',x,miu2,sigma1)';
            fprintf('Lognormal:\tmiu=%f\tlambda=%f\n',miu2,sigma1)
        case 'Rician'
            [s,sigma2] = estimateRicianParameters(x,pEmp,measure);
            pFit(:,i) = cdf('Rician',x,s,sigma2)';
        case 'Burr'
            [alpha,c1,k2] = estimateBurrParameters(x,pEmp,measure);
            pFit(:,i) = cdf('Burr',x,alpha,c1,k2)';
        case 'BirnbaumSaunders'
            [beta,gamma] = estimateBirnbaumSaundersParameters(x,pEmp,measure);
            pFit(:,i) = cdf('BirnbaumSaunders',x,beta,gamma)';
            fprintf('BirnbaumSaunders:\tbeta=%f\tgamma=%f\n',beta,gamma)
        case 'Rayleigh'
            [b1] = estimateRayleighParameters(x,pEmp,measure);
            pFit(:,i) = cdf('Rayleigh',x,b1)';
        case 'Weibull'
            [a,b2] = estimateWeibullParameters(x,pEmp,measure);
            pFit(:,i) = cdf('Weibull',x,a,b2)';
        case 'Noncentral Chi-square'
            [v,delta1] = estimateNCChi2Parameters(x,pEmp,measure);
            pFit(:,i) = cdf('Noncentral Chi-square',x,v,delta1)';
        case 'Noncentral F'
            [v11,v21,delta2] = estimateNCFParameters(x,pEmp,measure);
            pFit(:,i) = cdf('Noncentral F',x,v11,v21,delta2)';
        case  'F'
            [v12,v22] = estimateFParameters(x,pEmp,measure);
            pFit(:,i) = cdf('F',x,v12,v22)';
        case  'Nakagami'
            [miu3,omega] = estimateNakagamiParameters(x,pEmp,measure);
            pFit(:,i) = cdf('Nakagami',x,miu3,omega)';
        case  'Logistic'
            [miu4,sigma3] = estimateLogisticParameters(x,pEmp,measure);
            pFit(:,i) = cdf('Logistic',x,miu4,sigma3)';
        case  'LogLogistic'
            [miu5,sigma4] = estimateLogLogisticParameters(x,pEmp,measure);
            pFit(:,i) = cdf('LogLogistic',x,miu5,sigma4)';
        case  'Levy'
            %            [miu6,c2] = estimateLevyParameters(x,pEmp,0);
            [c2]=estimateLevyParameters0(x,pEmp,measure);
            pFit(:,i) = CDFOfLevy(x,0,c2);
    end
    pointCount = maxPointPloted;
    jump=ceil(size(pEmp,1)/pointCount);
    plot(pEmp(1:jump:size(pEmp,1)),pFit(1:jump:size(pEmp,1),i),legendList{i},'MarkerSize',10);
    hold on;
end
hold off;
xlabel('Empirical Probabilities');
ylabel('Fitted Probabilities');
legend([{'1:1 Line'},names], 'location','southeast');

% plot(x, pEmp,'r+')
% hold on
% plot(x, cdf('InverseGaussian',x,miu1,lambda))
% hold on
% legend([{'Experimental'},'inverse Gaussian'], 'location','southeast');
% xlabel('Pages');
% ylabel('Probility');
end

function [k, theta] = estimateGammaParameters(x,pEmp,measure)
a1=2;b1=1;
% wgt = 1 ./ sqrt(pEmp.*(1-pEmp));
obj = @(params) measure(gamcdf(x,exp(params(1)),exp(params(2)))-pEmp);
paramHat = fminsearch(obj,[log(a1),log(b1)]);

% options = optimoptions(@lsqnonlin,'Algorithm','Levenberg-Marquardt');
% paramHat = lsqnonlin(obj,[log(a1),log(b1)],[],[],options);

p = exp(paramHat);
k = p(1);
theta = p(2);
end



function [a1, b1] = estimateBetaParameters(x,pEmp,measure)
a1=2;b1=2;
% wgt = 1 ./ sqrt(pEmp.*(1-pEmp));
obj = @(params) measure(cdf('beta', x,exp(params(1)),exp(params(2)))-pEmp);
paramHat = fminsearch(obj,[log(a1),log(b1)]);

% options = optimoptions(@lsqnonlin,'Algorithm','Levenberg-Marquardt');
% paramHat = lsqnonlin(obj,[log(a1),log(b1)],[],[],options);

p = exp(paramHat);
a1 = p(1);
b1 = p(2);
end


function [miu, lambda] = estimateIGParameters(x,pEmp,measure)
a1=30;b1=17;
% wgt = 1 ./ sqrt(pEmp.*(1-pEmp));
obj = @(params) measure(cdf('InverseGaussian',x,exp(params(1)),exp(params(2)))-pEmp);
paramHat = fminsearch(obj,[log(a1),log(b1)]);

% options = optimoptions(@lsqnonlin,'Algorithm','Levenberg-Marquardt');
% paramHat = lsqnonlin(obj,[log(a1),log(b1)],[],[],options);


p = exp(paramHat);
miu = p(1);
lambda = p(2);
end

function [v, sigma] = estimateRicianParameters(x,pEmp,measure)
a1=30;b1=17;
% wgt = 1 ./ sqrt(pEmp.*(1-pEmp));
obj = @(params) measure(cdf('Rician',x,exp(params(1)),exp(params(2)))-pEmp);
paramHat = fminsearch(obj,[log(a1),log(b1)]);
p = exp(paramHat);
v = p(1);
sigma = p(2);
end


function [miu, sigma] = estimateLognParameters(x,pEmp,measure)
a1=30;b1=17;
% wgt = 1 ./ sqrt(pEmp.*(1-pEmp));
obj = @(params) measure(cdf('Lognormal',x,exp(params(1)),exp(params(2)))-pEmp);
paramHat = fminsearch(obj,[log(a1),log(b1)]);
p = exp(paramHat);
miu = p(1);
sigma = p(2);
end

function [beta,gamma] = estimateBirnbaumSaundersParameters(x,pEmp,measure)
a1=30;b1=17;
% wgt = 1 ./ sqrt(pEmp.*(1-pEmp));
obj = @(params) measure(cdf('BirnbaumSaunders',x,exp(params(1)),exp(params(2)))-pEmp);
paramHat = fminsearch(obj,[log(a1),log(b1)]);
p = exp(paramHat);
beta = p(1);
gamma = p(2);
end

function [alpha,c,k] = estimateBurrParameters(x,pEmp,measure)
alpha =10;
c=10;
k=10;
% wgt = 1 ./ sqrt(pEmp.*(1-pEmp));
obj = @(params) measure(cdf('Burr',x,exp(params(1)),exp(params(2)),exp(params(3)))-pEmp);
paramHat = fminsearch(obj,[log(alpha),log(c),log(k)]);
p = exp(paramHat);
alpha = p(1);
c = p(2);
k = p(3);
end

function [a,b] = estimateWeibullParameters(x,pEmp,measure)
a = 10;
b = 10;
% wgt = 1 ./ sqrt(pEmp.*(1-pEmp));
obj = @(params) measure(cdf('Weibull',x,exp(params(1)),exp(params(2)))-pEmp);
paramHat = fminsearch(obj,[log(a),log(b)]);
p = exp(paramHat);
a = p(1);
b = p(2);
end

function [v,delta] = estimateNCChi2Parameters(x,pEmp,measure)
v = 10;
delta = 10;
% wgt = 1 ./ sqrt(pEmp.*(1-pEmp));
obj = @(params) measure(cdf('Noncentral Chi-square',x,exp(params(1)),exp(params(2)))-pEmp);
paramHat = fminsearch(obj,[log(v),log(delta)]);
p = exp(paramHat);
v = p(1);
delta = p(2);
end

function [v1,v2,delta] = estimateNCFParameters(x,pEmp,measure)
v1=19;v2=19;delta=3;
% wgt = 1 ./ sqrt(pEmp.*(1-pEmp));
obj = @(params) measure(cdf('Noncentral F',x,exp(params(1)),exp(params(2)),exp(params(3)))-pEmp);
paramHat = fminsearch(obj,[log(v1),log(v2),log(delta)]);
p = exp(paramHat);
v1 = p(1);
v2 = p(2);
delta = p(3);
end

function [v1,v2] = estimateFParameters(x,pEmp,measure)
v1=19;v2=19;
% wgt = 1 ./ sqrt(pEmp.*(1-pEmp));
obj = @(params) measure(cdf('F',x,exp( params(1)),exp(params(2)))-pEmp);
paramHat = fminsearch(obj,[log(v1),log(v2)]);
p = exp(paramHat);
v1 = p(1);
v2 = p(2);
end

function [miu,omega] = estimateNakagamiParameters(x,pEmp,measure)
miu=5;omega=5;
% wgt = 1 ./ sqrt(pEmp.*(1-pEmp));
obj = @(params) measure(cdf('Nakagami',x,exp(params(1)),exp(params(2)))-pEmp);
paramHat = fminsearch(obj,[log(miu-0.5),log(omega)]);
p = exp(paramHat);
miu = p(1)+0.5;
omega = p(2);
end

function [miu,sigma] = estimateLogisticParameters(x,pEmp,measure)
miu=19;sigma=19;
%  wgt = 1 ./ sqrt(pEmp.*(1-pEmp));
obj = @(params) measure(cdf('Logistic',x,exp(params(1)),exp(params(2)))-pEmp);
paramHat = fminsearch(obj,[log(miu),log(sigma)]);
p = exp(paramHat);
miu = p(1);
sigma = p(2);
end


function [miu,sigma] = estimateLogLogisticParameters(x,pEmp,measure)
miu=19;sigma=19;
% wgt = 1 ./ sqrt(pEmp.*(1-pEmp));
obj = @(params) measure(cdf('LogLogistic',x,exp(params(1)),exp(params(2)))-pEmp);
paramHat = fminsearch(obj,[log(miu),log(sigma)]);
p = exp(paramHat);
miu = p(1);
sigma = p(2);
end

function [b] = estimateRayleighParameters(x,pEmp,measure)
b=19;
% wgt = 1 ./ sqrt(pEmp.*(1-pEmp));
obj = @(params) measure(cdf('Rayleigh',x,exp(params(1)))-pEmp);
paramHat = fminsearch(obj,log(b));
p = exp(paramHat);
b = p(1);
end

function [c] = estimateLevyParameters0(x,pEmp,measure)
c=6;
%  wgt = 1 ./ sqrt(pEmp.*(1-pEmp));
obj = @(params) measure(CDFOfLevy(x,0,exp(params(1)))-pEmp);
paramHat = fminsearch(obj,log(c));
p = exp(paramHat);
c = p(1);
end

function [miu sigma] = estimateExtremeValueParameters(x,pEmp,measure)
miu = 0;
sigma = 1;
obj = @(params) measure(evcdf(x,params(1),exp(params(2)))-pEmp);
paramHat = fminsearch(obj,[miu, sigma]);
miu = paramHat(1);
sigma = exp(paramHat(2));
end


function [miu,c] = estimateLevyParameters(x,pEmp,measure,threshold)
miu=-10;c=6;
if nargin == 3
    threshold = min(x)-0.0000001;
end
% wgt = 1 ./ sqrt(pEmp.*(1-pEmp));
obj = @(params) measure(CDFOfLevy(x,threshold/(1+exp(-params(1))),exp(params(2)))-pEmp);

paramHat = fminsearch(obj,[miu,log(c)]);
p = exp(paramHat);
miu=threshold/(1+exp(-paramHat(1)));
c = p(2);
end

function [p] = CDFOfLevy(x,miu, c)
p=x;
for i=1:size(p)
    p(i) = erfc(sqrt(c/(2*(x(i)-miu))));
end
end
