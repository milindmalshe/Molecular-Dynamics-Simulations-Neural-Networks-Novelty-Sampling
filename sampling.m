atoms=5

dat=load('data_UNsort_dihedralSign');
rTh_fConfig=dat.rTh_fConfig; 
clear dat;

temp=size(rTh_fConfig);
dimension=temp(1);
total=temp(2);

r_in=rTh_fConfig(1:(3*atoms-6),:);
f_in=rTh_fConfig((3*atoms-6)+1:dimension,:);
% clear data;

nn = load('NN-75-9_tansig-purelin&minmax');
NN5 = nn.net;
minr5 = nn.minr;
maxr5 = nn.maxr;
minf5 = nn.minf;
maxf5 = nn.maxf;

[rTh_fn,minrTh_f,maxrTh_f]=premnmx(rTh_fConfig);
rtr2=rTh_fn;
% r_tr=data.r_train;
% f_tr=data.f_train;
% val=data.val;
% % test=data.test;
% r_in=[r_tr val.P ];
% f_in=[f_tr val.T ];
% 
% 
% rtr2=  [r_tr;f_tr];
% rval2= [val.P;val.T];

%FOLLOWING IS TEMPORARY CHECK
%following config corresponds to experimental equilibrium config of 5 Si atoms
% rtr_check= [2.35166946731614 ; 2.35166946731614 ; 2.35166946731614 ; 2.35166946731614 ; 109.471220634491 ; 109.471220634491 ; 109.471220634491 ; 120 ; 120 ; -13.3965511648561];
%following config corresponds to G98 MP4(SDQ) 6-31G** energy minimum of 5 Si atoms
% rtr_check= [2.221 ; 2.221 ; 2.221 ; 2.221 ; 109.4712 ; 109.4712 ; 109.4712 ; 120 ; 120 ; -13.89782673];

% rtr2 = [rtr2 rtr_check];
%ABOVE IS TEMPORARY CHECK

temp=size(rtr2);
total=temp(2);

% new=rval2;
new=rtr2;
tempNew=size(new);
totalNew=tempNew(2);

countTotal=total;
% countTotal=300;

for i=24:countTotal 
    %     min(:,i) = new(:,i)-rtr2(:,1);
    if(i < countTotal)
        minDist_vector = rtr2(:,i)-rtr2(:,i+1);
		pos(i) = i+1;
        %         minDist(:,i) = sqrt(sum(minDist_vetcor(:,i).^2));
    else
        minDist_vector = rtr2(:,i)-rtr2(:,1);
		pos(i) = 1;
        %         minDist(:,i) = sqrt(sum(minDist_vetcor(:,i).^2));
    end
    minDist(i) = sqrt(sum(minDist_vector.^2));    
end



for i=24:countTotal
    for j=1:countTotal
        %         if((j ~= i) & (j ~= i+1))
        if(j ~= i)
            %             if(i < countTotal)
            minDist_vector = rtr2(:,i) - rtr2(:,j);
            %             else
            %                 minDist_vector = rtr2(:,i) - rtr2(:,1);
            %             end
            
            newDist = sqrt(sum(minDist_vector.^2));    
            
            if(newDist < minDist(i))
                minDist(i) = newDist;
                pos(i)=j;
            end
            
        end
    end
    i
end

p=1;% dimensionality =1, since dist is just one number
sigma=0.01;

% m=1;
% for i=0:0.01:100
%     xi(m)=i;
%     
%     sumj=0;
%     for j=1:countTotal
%         
%         sumj= sumj + exp(-1*(xi(m)-minDist(j)).^2/(2*sigma^2));
%     end
%     
%     fxi(m)= (1/((2*pi)^(p/2)*sigma^p)/countTotal)*sumj;
%     
%     m=m+1
% end

factor = 1/(sqrt(2*pi)*sigma*countTotal);

m=1;
for i=0:0.05:80
    fxi(m)= factor * sum(exp(-((i-minDist).^2)./(2*sigma^2)));
    m=m+1;
end

       