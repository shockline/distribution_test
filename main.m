data = load('data/imageDataAllSub.mat', 'imageDataAllSub');
data = data.imageDataAllSub;


[pValue,pFit,pEmp] = fitnessTest(data,'L1Sum',60);

