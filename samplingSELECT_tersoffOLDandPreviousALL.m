
atoms=5;

dat=load('data_UNsort_dihedralsign_tersoffOLD')
rTh_fConfig_tersoffOLD=dat.rThConfig_sort; 
temp=size(rTh_fConfig_tersoffOLD);
dimension=temp(1);
countTotalOLD=temp(2);
clear dat;

dat=load('data_UNsort_dihedralSign')
rTh_fConfig_ALL=dat.rTh_fConfig(1:(3*atoms-6),:);
temp=size(rTh_fConfig_ALL);
dimension=temp(1);
countTotalALL=temp(2);

[rTh_fn,minrTh_f,maxrTh_f]= premnmx(rTh_fConfig_ALL);
rTh_fn2= tramnmx(rTh_fConfig_tersoffOLD,minrTh_f,maxrTh_f)

for i=1:countTotalOLD 
    %     min(:,i) = new(:,i)-rtr2(:,1);
%     if(i < countTotal)
        minDist_vector = rTh_fn2(:,i) - rTh_fn(:,1);
		pos(i) = 1;
        %         minDist(:,i) = sqrt(sum(minDist_vetcor(:,i).^2));
%     else
%         minDist_vector = rtr2(:,i)-rtr2(:,1);
% 		pos(i) = 1;
%         %         minDist(:,i) = sqrt(sum(minDist_vetcor(:,i).^2));
%     end
    minDist(i) = sqrt(sum(minDist_vector.^2));    
    i
end



for i=1:countTotalOLD
    for j=1:countTotalALL
        %         if((j ~= i) & (j ~= i+1))
%         if(j ~= i)
            %             if(i < countTotal)
            minDist_vector = rTh_fn2(:,i) - rTh_fn(:,j);
            %             else
            %                 minDist_vector = rtr2(:,i) - rtr2(:,1);
            %             end
            
            newDist = sqrt(sum(minDist_vector.^2));    
            
            if(newDist < minDist(i))
                minDist(i) = newDist;
                pos(i)=j;
            end
            
%         end
    end
    i
end

minDist;
