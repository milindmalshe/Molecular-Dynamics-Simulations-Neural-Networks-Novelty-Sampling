

function[]= md()
clear;% ****THIS MODIFICATION IS MADE WHILE STORING THE CONFIGURATIONS 
mass = 27.9769; %atomic weight of Silicon, required to calculate acceleration from force: F = mass*accelrn 
temperature = 300; % temperature in KELVIN
% Boltzmann = 8.617547856e-5;
Boltzmann = 8.617402194e-5; %eV/K

% global delT
% delT = 0.05;


% global coord ;
%global x y z;

% global total numMov numPeriph numBound;
movAtom=0; boundAtom=0; periphAtom=0; surfaceAtom=0;

% global rCutoff;

clust_size=3;
type=[2,2,2];



totalITERATION=300;


coord=[ 0.000000        0.000000        0.000000;
	   -1.357737       1.357737        1.357737;
        1.357737        -1.357737       1.357737];

% coordtemp = load('coord_checkBackint');
% coord=coordtemp.coord;
   
x=coord(:,1); y=coord(:,2); z=coord(:,3);

 
   
   momentumBoltzmann = sqrt(2*mass*Boltzmann*temperature)/mass/sqrt(2);
%    momentumBoltzmann = sqrt(2*mass*Boltzmann*temperature)/mass/sqrt(3);
%    list=bondList(coord,numMov,total,rCutoff,movAtom); 
   rand('state',sum(100*clock));
   vel=zeros(total,3);
   for i=1:numMov
	   iMov=movAtom(i);
	   randNum = rand;
	   if(randNum < 0.5)
		   vel(iMov,1) = -1 * momentumBoltzmann;
	   else
		   vel(iMov,1) = momentumBoltzmann;
	   end
	   
	   randNum = rand;
	   if(randNum < 0.5)
		   vel(iMov,2) = -1 * momentumBoltzmann;
	   else
		   vel(iMov,2) = momentumBoltzmann;
	   end
	   
	   randNum = rand;
	   if(randNum < 0.5)
		   vel(iMov,3) = -1 * momentumBoltzmann;
	   else
		   vel(iMov,3) = momentumBoltzmann;
	   end
   end
   
   
   
   
%    vel=zeros(total,3);
% vel=zeros(numMov,3);
   accelrn=zeros(total,3);
% accelrn=zeros(numMov,3);
 
% [tersoff_PE,force] = tersoffSi3_PeriodicBoundaryCond(coord,total,numMov,numPeriph,numBound,numSurface,movAtom,periphAtom,boundAtom,surfaceAtom);
[tersoff_PE,force] = NNG98_2_coordStore_PBC_4ForceComponentOnALLatomsInList(coord,total,numMov,numPeriph,numBound,numSurface,movAtom,periphAtom,boundAtom,surfaceAtom);
% [tersoff_PE,force] = NNG98_2_coordStore(coord,total,numMov,numPeriph,numBound,movAtom,periphAtom,boundAtom);

accelrn = -1.*force./mass;

% KE = zeros(total,1);
KE=0;
for i=1:total
	velTot(i)= sqrt(sum(vel(i,:).^2));
end

for i=1:total
	KE= KE + 1/2*mass.*(velTot(i).^2);
end
% hndl = scatter3(coord(:,1),coord(:,2),coord(:,3));
% set(hndl, 'EraseMode', 'xor', 'MarkerSize',10)

for count=1:totalITERATION
		
%     if (rem(count,5)==0)
%         coord(1,1)  = coord(1,1)  + 0.1*rand;   coord(1,2) = coord(1,2)  + 0.1*rand;  coord(1,3) =  coord(1,3)  + 0.1*rand;
%         coord(71,1) = coord(71,1) + 0.1*rand;  coord(71,2) = coord(71,2) + 0.1*rand;  coord(71,3) = coord(71,3) + 0.1*rand;
%         coord(73,1) = coord(73,1) + 0.1*rand;  coord(73,2) = coord(73,2) + 0.1*rand;  coord(73,3) = coord(73,3) + 0.1*rand;
%         coord(82,1) = coord(82,1) + 0.1*rand;  coord(82,2) = coord(82,2) + 0.1*rand;  coord(82,3) = coord(82,3) + 0.1*rand;
%         coord(92,1) = coord(92,1) + 0.1*rand;  coord(92,2) = coord(92,2) + 0.1*rand;  coord(92,3) = coord(92,3) + 0.1*rand;
%         
% %         coord(14,2)= coord(14,2)+ 0.1*rand;
% %         coord(15,3)= coord(15,3)+ 0.1*rand;
% %         coord(16,1)= coord(16,1)+ 0.1*rand;
% %         coord(17,2)= coord(17,2)+ 0.1*rand;
%         
%         [coord,vel,accelrn,tersoff_PE,KE] = velverlet(count,coord,total,numMov,numPeriph,numBound,movAtom,periphAtom,boundAtom,vel,accelrn,mass);
        	 
        coord;
%     else
        [coord,vel,accelrn,tersoff_PE,KE] = velverlet_compPhys(count,coord,total,numMov,numPeriph,numBound,numSurface,movAtom,periphAtom,boundAtom,surfaceAtom,vel,accelrn,mass);
%     end
    count
        
	if(count ==10)
		coord;
	end
% 	set(hndl,'XData',coord(:,1),'YData',coord(:,2),'ZData',coord(:,3))
% 	drawnow
% 	kineticE=sqrt(sum(KE.^2));
% tersoff_PE;
    totalE(count) = KE+tersoff_PE;
    totalE(count)	
% 	if(rem(count,100) == 0)
% 		tersoff_PE
% 		KE
% 		totalE= KE + tersoff_PE
% 		count
% 	end
%FOLLOWING CODE IS TO WRITE COORDINATES FOR ANIMATION FILE
if(rem(count,300) == 0)
        fileCount=num2str(count/300)
        if(str2num(fileCount) < 10)
           fileWRKR = strcat('wrkr037.f0',fileCount);
           fileTOLR = strcat('tolr037.f0',fileCount);
        else
           fileWRKR = strcat('wrkr037.f',fileCount);
           fileTOLR = strcat('tolr037.f',fileCount);
        end

        fidWRKR=fopen(fileWRKR,'w');
        fidTOLR=fopen(fileTOLR,'w');

        for (i=1:numMov)
                fprintf(fidWRKR,'%f\t%f\t%f',coord(movAtom(i),1),coord(movAtom(i),2),coord(movAtom(i),3));
                fprintf(fidWRKR,'\n');
        end

        fclose(fidWRKR);
        fclose(fidTOLR);
 end
%ABOVE CODE IS TO WRITE COORDINATES FOR ANIMATION FILE

% %FOLLOWING CODE IS TO RESET THE VELOCITY OF PERIPHERAL ATOMS
% for i=1:numPeriph
%     iPeriph=periphAtom(i);
%     randNum = rand;
%     if(randNum < 0.5)
%         vel(iPeriph,1) = sqrt(1-0.1047)*vel(iPeriph,1) + sqrt(0.1047) * momentumBoltzmann;
%     else
%         vel(iPeriph,1) = sqrt(1-0.1047)*vel(iPeriph,1) - sqrt(0.1047) * momentumBoltzmann;
%     end
%     
%     randNum = rand;
%     if(randNum < 0.5)
%         vel(iPeriph,2) = sqrt(1-0.1047)*vel(iPeriph,2) + sqrt(0.1047) * momentumBoltzmann;
%     else
%         vel(iPeriph,2) = sqrt(1-0.1047)*vel(iPeriph,2) - sqrt(0.1047) * momentumBoltzmann;
%     end
%     
%     randNum = rand;
%     if(randNum < 0.5)
%         vel(iPeriph,3) = sqrt(1-0.1047)*vel(iPeriph,3) + sqrt(0.1047) * momentumBoltzmann;
%     else
%         vel(iPeriph,3) = sqrt(1-0.1047)*vel(iPeriph,3) - sqrt(0.1047) * momentumBoltzmann;
%     end
% end
% %ABOVE CODE IS TO RESET THE VELOCITY OF PERIPHERAL ATOMS
end

plot(1:500,totalE,'.')
  
