x = 1 : 800;
miu = 30.645818;
lambda = 17.113267;
y  = (cdf('InverseGaussian',x,miu,lambda)-cdf('InverseGaussian',x-1,miu,lambda))./(1-cdf('InverseGaussian',x-1,miu,lambda));
plot(x,1-y);
hold on;

z = imageDataAll(:,1:2);
zx = z(:,1);
zy = z(:,2);
emp = cumsum(zy)/sum(zy);
zy = (emp(2:end)-emp(1:end-1))./(1-emp(1:end-1));
plot(zx(2:end),1-zy);
hold off;