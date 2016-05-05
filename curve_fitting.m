method='L1Sum'
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
obj = @(params) measure(CDFOfLevy(x,0,exp(params(1)))-pEmp);
paramHat = fminsearch(obj,log(c));
p = exp(paramHat);
c = p(1);