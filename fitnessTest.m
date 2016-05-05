%% data
% 1st column: sample
% 2nd column: sample probility
% 3rd column: fruency
%% Measure
% switch method
%     case 'L1Sum'
%         measure = @(params) sum(abs(params));
%     case 'L1Max'
%         measure = @(params) max(abs(params));
%     case 'L2Sum'
%         measure = @(params) sum((params).^2);
%     otherwise
%         measure = @(params) sum((params).^2);
% end
function [pValue,pFit,pEmp] = fitnessTest(data,method,maxPointPloted)
if nargin == 1
    method = 'L2Sum';
end

if nargin == 2
    maxPointPloted = 80;
end

%% Work
data=sortrows(data, 1);
% if sum(data(:,2))==0
    data(:,2)=data(:,3)/sum(data(:,3));
% end

[pEmp, names, pFit]= CDFFit([data(:,1), cumsum(data(:,2))], method, maxPointPloted);
probCount = size(pFit,2);
pValue = zeros(1, probCount);

for i = 1:probCount
    pValue(1,i) = ksTest(pEmp,pFit(:,i));
end

% for i = 1:probCount
%     pValue(2,i) = adTest(pFit(1:24,i), data(1:24,1) , data(1:24,3) );
% end

%% Set KS-test confidence level
% \alpha	    0.10	0.05	0.025	0.01	0.005	0.001
% c({\alpha})	1.22	1.36	1.48	1.63	1.73	1.95
alpha = 0.05;
level = 1.36  *  sqrt(2/sum(data(:,3)));

%% Output
fprintf('Compare with: %f(alpha = %f)\n', level, alpha);
for i = 1:probCount
    fprintf('%20s\t%f',names{i},pValue(i));
    if pValue(i) > level
        fprintf('\tRejcted\n');
    else
        fprintf('\tAccepted\n');
    end
end

pValue = zeros(1, probCount);
hValue = zeros(1, probCount);
for i = 1:probCount
    [n,p]=size(pEmp);
    y = pEmp;
    y_hat=pFit(:,i);
    
    [muhat,sigmahat] = normfit(y-y_hat);
%     p=1;
     SSR   =     (y_hat-mean(y))'*(y_hat-mean(y)); % 回归平方和
     SSE = (y-y_hat)'*(y-y_hat);
%      SSTO   =     (y-mean(y))'*(y-mean(y));         % 完全平方和

     statistic = SSR/(SSR+SSE);
    
%      F     =     (SSR/p)/(SSE/(n-1-p));            % F
%      sig   =     1-fcdf(F,p,n-p-1);
%     m = mean(pEmp);
%     p1 = pEmp -m;
%     p2 = pFit(:,i) -m;
%     RR = sqrt(sum(p1.^2)/sum(p2.^2))

     hValue(1,i) = statistic;
     pValue(1,i) = sigmahat;
%      [hValue(1,i) pValue(1,i)] = vartest2(pEmp,pFit(:,i),'Alpha',0.95);
end
%% F-test
%fprintf('Compare with: %f(alpha = %f)\n', level, alpha);
fprintf('F-test\n', level, alpha);
for i = 1:probCount
    fprintf('%20s\t%f\t%f',names{i},hValue(i),pValue(i));
    if hValue(i) == 1
        fprintf('\tRejcted\n');
    else
        fprintf('\tAccepted\n');
    end
end

end
%% KS-test
function [p] = ksTest(realCumP, fitCumP)
count=size(realCumP,1);
p1 = max(abs(realCumP-fitCumP));
p2 = max(abs(realCumP(1:count-1)-fitCumP(2:count)));
p=max(p1,p2);
end

%% AD-test
%{
function [p] = adTest(fitCumP,x,numberofx)
    n = sum(numberofx);
    cumValue = cumsum(numberofx);
    S = 0;
    lowerbound = 1;
    for line=1:n
        if(line > cumValue(lowerbound))
            lowerbound = lowerbound + 1;
        end
        y0 = fitCumP(lowerbound);
        S = S +  (2*line-1)*log(y0) + (2*(n-line)+1)*log (1-y0+0.0000001)   ;
    end
    AD = -n - S/n;
    p=AD;
    end
%}