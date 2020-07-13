function fFitness=calculateFitness(fObjV)
fFitness=zeros(size(fObjV));
ind=find(fObjV>=0);
fFitness(ind)=1./(fObjV(ind)*0.5+1);   % 设为0.5效果是最好的，对于布谷鸟来说
% fFitness(ind) = exp(-fObjV(ind)*0.05);
ind=find(fObjV<0);  % 怎么会有小于0的呢？!
fFitness(ind)=1+abs(fObjV(ind));
