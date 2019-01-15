
d=load('minDist_CONCATANATE_sqrt_NEW_DEL_RECALCULATE')
% hist(d.minDist_CONCATANATE_sqrt_NEW_DEL_RECALCULATE,50)
tempSize=size(d.minDist_CONCATANATE_sqrt_NEW_DEL_RECALCULATE);
total=tempSize(2);

p=1;
sigma=0.001;
factor = 1/(sqrt(2*pi)*sigma*(total));
m=1;
for i=0:0.01:0.2
    fxi(m)= factor * sum(exp(-((i-d.minDist_CONCATANATE_sqrt_NEW_DEL_RECALCULATE).^2)./(2*sigma^2)));
    m=m+1;
end
% plot(0:.01:.2,fxi)
dNEW=load('minDist_FINAL')
tempSize=size(dNEW.minDistNEW);
totalNEW=tempSize(2);

factorNEW = 1/(sqrt(2*pi)*sigma*(totalNEW));
m=1;

T1=0.1;

accept1=0; accept2=0; reject=0;

for i=1:totalNEW
    fxiNEW= factor * sum(exp(-((dNEW.minDistNEW(i)-d.minDist_CONCATANATE_sqrt_NEW_DEL_RECALCULATE).^2)./(2*sigma^2)));
    
    if(fxiNEW <= T1*max(fxi) & fxiNEW ~= 0 & dNEW.minDistNEW(i) > 0.07)
        accept1 = accept1+1;
    elseif(fxiNEW./max(fxi) <= rand*rand*(1-T1) & fxiNEW ~= 0 & dNEW.minDistNEW(i) > 0.07)
        accept2=accept2+1;
    else
        reject=reject+1;
    end
end

T1=0.1;

accept1=0; accept2=0; reject=0;

for i=1:totalNEW
    if(dNEW.minDistNEW(i) <= T1*max(fxi) & dNEW.minDistNEW(i) ~= 0)
        accept1 = accept1+1;
    elseif(dNEW.minDistNEW(i)./max(fxi) <= rand*(1-T1) & dNEW.minDistNEW(i) ~= 0)
        accept2=accept2+1;
    else
        reject=reject+1;
    end
end


% for i=1:totalNEW
%     if(dNEW.minDistNEW(i) <= T1*max(fxi) & dNEW.minDistNEW(i) ~= 0)
%         accept1 = accept1+1;
%     elseif(dNEW.minDistNEW(i)./max(fxi) <= rand*(1-T1) & dNEW.minDistNEW(i) ~= 0)
%         accept2=accept2+1;
%     else
%         reject=reject+1;
%     end
% end
