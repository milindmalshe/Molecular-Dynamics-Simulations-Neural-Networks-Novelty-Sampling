atoms=5

dat=load('data_CONCATANATE_sqrt_NEW_DEL');
rTh_fConfig=dat.data_CONCATANATE_sqrt_NEW_DEL; 
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

dat=load('minDist_CONCATANATE_sqrt_NEW_DEL');
minDist_sqrt=dat.minDist_CONCATANATE_sqrt_NEW_DEL;

% for i=24:countTotal 
%     %     min(:,i) = new(:,i)-rtr2(:,1);
%     if(i < countTotal)
%         minDist_vector = rtr2(:,i)-rtr2(:,i+1);
% 		pos(i) = i+1;
%         %         minDist(:,i) = sqrt(sum(minDist_vetcor(:,i).^2));
%     else
%         minDist_vector = rtr2(:,i)-rtr2(:,1);
% 		pos(i) = 1;
%         %         minDist(:,i) = sqrt(sum(minDist_vetcor(:,i).^2));
%     end
%     minDist(i) = sqrt(sum(minDist_vector.^2));    
% end
% 
% 
% 
% for i=24:countTotal
%     for j=1:countTotal
%         %         if((j ~= i) & (j ~= i+1))
%         if(j ~= i)
%             %             if(i < countTotal)
%             minDist_vector = rtr2(:,i) - rtr2(:,j);
%             %             else
%             %                 minDist_vector = rtr2(:,i) - rtr2(:,1);
%             %             end
%             
%             newDist = sqrt(sum(minDist_vector.^2));    
%             
%             if(newDist < minDist(i))
%                 minDist(i) = newDist;
%                 pos(i)=j;
%             end
%             
%         end
%     end
%     i
% end
% 
% p=1;% dimensionality =1, since dist is just one number
% sigma=0.01;
% 
% % m=1;
% % for i=0:0.01:100
% %     xi(m)=i;
% %     
% %     sumj=0;
% %     for j=1:countTotal
% %         
% %         sumj= sumj + exp(-1*(xi(m)-minDist(j)).^2/(2*sigma^2));
% %     end
% %     
% %     fxi(m)= (1/((2*pi)^(p/2)*sigma^p)/countTotal)*sumj;
% %     
% %     m=m+1
% % end
% 
% factor = 1/(sqrt(2*pi)*sigma*countTotal);
% 
% m=1;
% for i=0:0.05:80
%     fxi(m)= factor * sum(exp(-((i-minDist).^2)./(2*sigma^2)));
%     m=m+1;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SELECT SAMPLING%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

datNEW=load('data_ALL_DEL');
rTh_fConfigNEW=datNEW.data_ALL_DEL; 
clear datNEW;

temp=size(rTh_fConfigNEW);
dimensionNEW=temp(1);
totalNEW=temp(2);

% r_inNEW=rTh_fConfigNEW(1:(3*atoms-6),total+1:totalNEW);
% f_inNEW=rTh_fConfigNEW((3*atoms-6)+1:dimensionNEW,total+1:totalNEW);
% 
% f_inNEW_sim=postmnmx(sim(NN5,tramnmx(r_inNEW,minr5,maxr5)),minf5,maxf5);
% 
% rTh_fConfigNEW_sim = [r_inNEW;f_inNEW_sim];
% 
% temp_sim=minmax(rTh_fConfigNEW_sim);
% minrTh_f_sim=temp_sim(:,1);
% maxrTh_f_sim=temp_sim(:,2);

[rTh_fnNEW]=tramnmx(rTh_fConfigNEW,minrTh_f(1:9),maxrTh_f(1:9));
rtr2NEW=rTh_fnNEW;

minDistNEW=zeros(1,totalNEW);
m=0
for i=total+1:totalNEW
    %     min(:,i) = new(:,i)-rtr2(:,1);
	%     if(i < countTotal)
	% 		%         minDist_vector = rtr2(:,i)-rtr2(:,i+1);
	% 		m=m+1;
	% 		minDist_vectorNEW = rtr2NEW(:,i)-rtr2(:,1);
	% 		% 		pos(i) = i+1;
	% 		pos(i) = 1;
	%         %         minDist(:,i) = sqrt(sum(minDist_vetcor(:,i).^2));
	%     else
		m=m+1;
        minDist_vectorNEW = rtr2NEW(:,m)-rtr2(1:9,1);
		pos(m) = 1;
        %         minDist(:,i) = sqrt(sum(minDist_vetcor(:,i).^2));
		%     end
    minDistNEW(m) = sqrt(sum(minDist_vectorNEW.^2));    
end


m=0;
for i=1:34371
	m=m+1;
    for j=1:total
        %         if((j ~= i) & (j ~= i+1))
        if(j ~= i)
            %             if(i < countTotal)
            minDist_vectorNEW = rtr2NEW(:,m) - rtr2(1:9,j);
            %             else
            %                 minDist_vector = rtr2(:,i) - rtr2(:,1);
            %             end
            
            newDistNEW = sqrt(sum(minDist_vectorNEW.^2));    
            
            if(newDistNEW < minDistNEW(m))
                minDistNEW(m) = newDistNEW;
                pos(m)=j;
            end
            
        end
    end
	%     i
	m
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
minDistNEW=load('minDistNEW');
minDistNEW=minDistNEW.minDistNEW;

factor = 1/(sqrt(2*pi)*sigma*(total+totalNEW));

m=1;
for i=0:0.01:0.25
    fxi(m)= factor * sum(exp(-((i-minDistNEW).^2)./(2*sigma^2)));
    m=m+1;
end
