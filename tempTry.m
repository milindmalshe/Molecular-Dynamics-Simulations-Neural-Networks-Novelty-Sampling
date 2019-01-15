%This program reads x,y,z configurations from separate x,y,z files and writes multiple input files for Gaussian 98 with multiple configurations separated by '--Link1--' command in G98
%The input file G98 first performs Hartree Fock calculations and then MP4 calculations, TO AVOIND CONVERGENCE PROBLEM

% total=19683;%Total number of configurations (cells) in the file
atoms=5;

x=xlsread('x5_coordSurface_Tersoff_10thselect.xls');%read x coordinates
y=xlsread('y5_coordSurface_Tersoff_10thselect.xls');%read y coordinates
z=xlsread('z5_coordSurface_Tersoff_10thselect.xls');%read z coordinates
x(:,1)=1.5.*x(:,1); x(:,2)=1.5.*x(:,2); y(:,1)=1.5.*y(:,1); y(:,2)=1.5.*y(:,2); z(:,1)=1.5.*z(:,1); z(:,2)=1.5.*z(:,2);
temp = size(x);
total = temp(1);%Total number of configurations (cells) in the file
total=1;
numConfig=100;% Number of configurations in a file
fileNumMax=total/numConfig;% calculate number of files to be generated = total number of configurations/Number of configurations per file 

r= zeros(atoms,atoms,total);
theta= zeros(atoms,atoms,atoms,total);
dihedralTheta= zeros(atoms,1,atoms,atoms,total);

for i=1:total
    for j=1:(atoms-1)
        r(1,j+1,i)=sqrt(x(i,j).^2+y(i,j).^2+z(i,j).^2);%calculate distances of other atoms from atom 1 
    end;
end;
    
order_sort=zeros(4,total);
for i=1:total
    order_sort(2,i)=2;
    order_sort(3,i)=3;
    order_sort(4,i)=4;
    order_sort(5,i)=5;
end

for i=1:total
    for j=1:(atoms)
        for k=2:(atoms)-j
            if(abs(r(1,k,i)-2.3517) > abs(r(1,k+1,i)-2.3517))
                temp=r(1,k,i);
                r(1,k,i)=r(1,k+1,i);
                r(1,k+1,i)=temp;
                
                temp_order=order_sort(k,i);
                order_sort(k,i)=order_sort(k+1,i);
                order_sort(k+1,i)=temp_order;
                
            end
        end
    end
end

x_NEW=zeros(total,4);
y_NEW=zeros(total,4);
z_NEW=zeros(total,4);
for i=1:total
	for j=2:atoms
		x_NEW(i,j-1)=x(i,order_sort(j,i)-1);
		y_NEW(i,j-1)=y(i,order_sort(j,i)-1);
		z_NEW(i,j-1)=z(i,order_sort(j,i)-1);
	end
end

x=x_NEW;
y=y_NEW;
z=z_NEW;

for i=1:total
    for j=1:(atoms-1);
        for k=1:(atoms-1);
            if(k>j)
%                 r(j,k,i)=(x_r(i,j)-x_r(i,k));
                                r(j+1,k+1,i)=sqrt((x(i,j)-x(i,k)).^2+(y(i,j)-y(i,k)).^2+(z(i,j)-z(i,k)).^2);%calculate distances of other atoms with each other
                
            end;
        end;
    end;
end;


for i=1:total
    for j=1:atoms
        for k=1:atoms
            if(k~=j & k<j)
                r(j,k,i)=r(k,j,i);%Make the r matrix symmetric
            end % end if
        end % end for k
    end % end for j
end % end for i

% for i=1:total
%     for j=1:atoms-1
%         [th(i,j),phi(i,j),r_(i,j)]=cart2sph(x(i,j),y(i,j),z(i,j));% convert cartesian coordinates to spherical coordinates using cart2sph, to get angle th and phi
%         th(i,j)=th(i,j).*180.0./pi;
%     end;
% end;


for i=1:total
    for j=1:atoms
        for k=1:atoms
            for m=1:atoms
                if (k~=j & m~=k & m~=j)
                    cosine(j,k,m,i)=(r(j,k,i).^2+r(k,m,i).^2-r(j,m,i).^2)./(2.*r(j,k,i).*r(k,m,i));% consine rule ** k-ATOM IS THE CENTRAL ATOM FOR COSINE RULE ** 
                    theta(j,k,m,i)=acos(cosine(j,k,m,i)).*180.0./pi;% calculate by taking acos and then cnovert the angle from radians to degrees
                end %end if
            end % end for m
        end %end for k
    end %end for j
end %end for i


for i=1:total
    for j=1:atoms
        for k=1:atoms
            for m=1:atoms
                if(j~=1 & k~=1 & m~=1 & m~=j & m~=k & k~=j)
                    x_n(1) = x(i,m-1) - x(i,j-1);
                    x_n(2) = x(i,j-1) - 0;    %ACTUALLY THIS IS x(2)-x(1) but since x(1)=0
                    x_n(3) = 0 - x(i,k-1);    %ACTUALLY THIS IS x(1)-x(3) but since x(1)=0
                    
                    y_n(1) = y(i,m-1) - y(i,j-1);
                    y_n(2) = y(i,j-1) - 0;    %ACTUALLY THIS IS y(2)-y(1) but since x(1)=0
                    y_n(3) = 0 - y(i,k-1);    %ACTUALLY THIS IS y(1)-y(3) but since x(1)=0
                    
                    z_n(1) = z(i,m-1) - z(i,j-1);
                    z_n(2) = z(i,j-1) - 0;    %ACTUALLY THIS IS z(2)-z(1) but since x(1)=0
                    z_n(3) = 0 - z(i,k-1);    %ACTUALLY THIS IS z(1)-z(3) but since x(1)=0
                    
                    %TO CALCULATE DIHEDRAL ANGLE, FIRST CALCULATE EQUATIONS OF NORMAL TO THE 2 PLANES FORMED viz. 3-1-2 AND 1-2-4
                    n_1x = y_n(2)*z_n(1) - y_n(1)*z_n(2);%calculate the X components of 1st normal vector to the plane 1-2-4
                    n_1y = z_n(2)*x_n(1) - x_n(2)*z_n(1);
                    n_1z = x_n(2)*y_n(1) - y_n(2)*x_n(1);
                    
                    n_2x = y_n(3)*z_n(2) - z_n(3)*y_n(2);%calculate the X components of 1st normal vector to the plane 3-1-2
                    n_2y = z_n(3)*x_n(2) - x_n(3)*z_n(2);
                    n_2z = x_n(3)*y_n(2) - y_n(3)*x_n(2);
                    
                    n_1 = sqrt(n_1x^2 + n_1y^2 + n_1z^2);%CALCULATE THE MAGNITUDE OF NORMAL VECTOR, GIVEN BY sqrt OF squares OF X,Y,Z COMPONENTS
                    n_2 = sqrt(n_2x^2 + n_2y^2 + n_2z^2);
                    
                    n_1DOTn_2 = n_1x*n_2x + n_1y*n_2y + n_1z*n_2z;%DOT PRODUCT OF THE 2 NORMAL VECTORS GIVEN BY SUM OF PRODUCTS OF X,Y,Z COMPONENTS OF 2 VECTORS
                    n_1_n_2 = n_1*n_2;%PRODUCT OF MAGNITUDES OF THE 2 NORMAL VETORS
                    
                    cosDihedral = n_1DOTn_2 / n_1_n_2; %cosine of dihedral angle
                    dihedralTheta(m,1,j,k,i) = acos(cosDihedral)*180/pi;%DIHEDRAL ANGLE
                end
            end
        end
    end
end


% 
% %FOLLOWING IS OLD CODE TO CALCULATE DIHEDRAL ANGLE
% for i=1:total% convert cartesian coordinate system to spherical coordiante system
%     for j=1:(atoms-1)
%         [th(i,j),phi(i,j),rsph(i,j)]=cart2sph(x(i,j),y(i,j),z(i,j));
%     end
% end
% 
% 
% % for i=1:total% ROTATE COORDINATE AXES SO THAT Z AXIS IS ALONG BOND '1-2', For EACH CONFIGURATION. 
% %              % THIS REQUIRES 1st ROTATION ABOUT Z-AXIS BY ANGLE (90-th), and 2nd ROTATION ABOUT X-AXIS BY (90-phi),
% %              % ....th & phi are calculated using cartesian to spherical transformation
% %     for j=1:(atoms-1)
% %         x_new2(i,j)=(x(i,j).*cos((pi./2)-th(i,1))-y(i,j).*sin((pi./2)-th(i,1)));
% %         y_new2(i,j)=(x(i,j).*sin((pi./2)-th(i,1))+y(i,j).*cos((pi./2)-th(i,1))).*cos((pi./2)-phi(i,1))-(z(i,j).*sin((pi./2)-phi(i,1)));
% %         z_new2(i,j)=(x(i,j).*sin((pi./2)-th(i,1))+y(i,j).*cos((pi./2)-th(i,1))).*sin((pi./2)-phi(i,1))-(z(i,j).*cos((pi./2)-phi(i,1)));
% %     end;
% % end;
% 
% 
% 
% 
% for i=1:total% ROTATE COORDINATE AXES SO THAT Z AXIS IS ALONG BOND '1-2', For EACH CONFIGURATION. 
%              % THIS REQUIRES 1st ROTATION ABOUT Z-AXIS BY ANGLE (90-th), and 2nd ROTATION ABOUT X-AXIS BY (90-phi),
%              % ....th & phi are calculated using cartesian to spherical transformation
%     for j=1:(atoms-1)
%         for k=1:(atoms-1)
% %             if(k ~= j)
%                 x_new(i,k,j)=(x(i,k).*cos((pi./2)-th(i,j))-y(i,k).*sin((pi./2)-th(i,j)));
%                 y_new(i,k,j)=(x(i,k).*sin((pi./2)-th(i,j))+y(i,k).*cos((pi./2)-th(i,j))).*cos((pi./2)-phi(i,j))-(z(i,k).*sin((pi./2)-phi(i,j)));
%                 z_new(i,k,j)=(x(i,k).*sin((pi./2)-th(i,j))+y(i,k).*cos((pi./2)-th(i,j))).*sin((pi./2)-phi(i,j))-(z(i,k).*cos((pi./2)-phi(i,j)));
% %             end
%         end
%     end
% end
%          
% 
% for i=1:total%CALCULATE XYZ COMPONENTS OF THE NORMAL TO PLANES 1-2-atom FOR EACH CONFIGURATION
%     for j=1:atoms
%         for k=1:atoms
%             if(j~=1 & k~=1 & k~=j)
%                 xComp(k,1,j,i)=     (y_new(i,k-1,j-1).*z_new(i,j-1,j-1))-(y_new(i,j-1,j-1).*z_new(i,k-1,j-1));
%                 yComp(k,1,j,i)=-1.*((x_new(i,k-1,j-1).*z_new(i,j-1,j-1))-(x_new(i,j-1,j-1).*z_new(i,k-1,j-1)));
%                 zComp(k,1,j,i)=     (x_new(i,k-1,j-1).*y_new(i,j-1,j-1))-(x_new(i,j-1,j-1).*y_new(i,k-1,j-1));
%                 
%                 
%                 
%                 
% %                 xComp2(k,1,2,i)=     (y_new(i,k-1).*z_new(i,2-1))-(y_new(i,2-1).*z_new(i,k-1));
% %                 yComp2(k,1,2,i)=-1.*((x_new(i,k-1).*z_new(i,2-1))-(x_new(i,2-1).*z_new(i,k-1)));
% %                 zComp2(k,1,2,i)=     (x_new(i,k-1).*y_new(i,2-1))-(x_new(i,2-1).*y_new(i,k-1));
% 
% 
%                 [th_12_Z(k,1,j,i),phi_12_Z(k,1,j,i),rsph_12_Z(k,1,j,i)]=cart2sph(xComp(k,1,j,i),yComp(k,1,j,i),zComp(k,1,j,i));%FIND th_12_Z in NEUMAN Projection
% 
%                 %             theta_12_Z(j,1,2,i)=atan(yComp(j,1,2,i)./xComp(j,1,2,i)); % theta_12_Z IS THE ANGLE BETWEEN ATOMS j-1-2 at 1, 
%                 % WHILE LOOKING ALONG BOND 1-2 BEING Z-AXIS
%                 
%                 %             if(theta_12_Z(j,1,2,i)<0.0)
%                 %                 theta_12_Z(j,1,2,i)=360-abs(theta_12_Z(j,1,2,i));% if angle is -ve (i.e anticlockwise as per Gaussian notation) the make +ve by taking angle=+ve (360-angle)
%                 %             end
%             end
%         end
%     end
% end
% 
% 
% for i=1:total% CALCULATE THE DIHEDRAL ANGLE
%     for j=1:atoms
%         for k=1:atoms
%             for m=1:atoms
%                 if(j~=1 & k~=1 & m~=1 & m~=j & m~=k & k~=j)
%                     
%                     numerator=xComp(k,1,j,i).*(xComp(m,1,j,i))+yComp(k,1,j,i).*(yComp(m,1,j,i))+zComp(k,1,j,i).*(zComp(m,1,j,i));
%                     denominator=(sqrt(xComp(k,1,j,i).^2+yComp(k,1,j,i).^2+zComp(k,1,j,i).^2)).*(sqrt(xComp(m,1,j,i).^2+yComp(m,1,j,i).^2+zComp(m,1,j,i).^2));
%                     
% %                     dihedralTheta(j,1,2,k,i)=acos((numerator./denominator)).*180.0./pi;
%                     dihedralTheta(k,1,j,m,i)=acos((numerator./denominator)).*180.0./pi;
%                     
%                     %                 if(i==1000)
%                     %                     temp=0
%                     %                 end
%                     
% %                     cosAlpha=xComp(j,1,2,i).*x_new(i,k-1)+yComp(j,1,2,i).*y_new(i,k-1)+zComp(j,1,2,i).*z_new(i,k-1);%cosAPLHA is > 1, SO DONOT EVALUATE ANGLE alpha
%                     cosAlpha=xComp(k,1,j,i).*x_new(i,m-1,j-1)+yComp(k,1,j,i).*y_new(i,m-1,j-1)+zComp(k,1,j,i).*z_new(i,m-1,j-1);%cosAPLHA is > 1, SO DONOT EVALUATE ANGLE alpha
%                                                                                                                     %CONSIDER ONLY THE SIGN OF cosAlpha TO DETERMINE
%                                                                                                                     %THE SIGN OF THE DIHEDRAL ANGLE
%                     
%                     %                 cosAlpha=xComp(j,1,2,i).*(x_new(i,k-1)-x_new(i,1))+yComp(j,1,2,i).*(y_new(i,k-1)-y_new(i,1))+zComp(j,1,2,i).*(z_new(i,k-1)-z_new(i,1));
%                     
%                     
%                     if(cosAlpha < 0.0)
%                         dihedralTheta(k,1,j,m,i)= (-1).*dihedralTheta(k,1,j,m,i);
%                     end
%                     
%                     
%                     
%                     
%                     
%                     xCompRotate_j=xComp(k,1,j,i).*cos((-1).*th_12_Z(k,1,j,i))-yComp(k,1,j,i).*sin(th_12_Z(k,1,j,i));
%                     yCompRotate_j=xComp(k,1,j,i).*sin((-1).*th_12_Z(k,1,j,i))+yComp(k,1,j,i).*cos(th_12_Z(k,1,j,i));
%                     
%                     xCompRotate_k=xComp(m,1,j,i).*cos((-1).*th_12_Z(k,1,j,i))-yComp(m,1,j,i).*sin(th_12_Z(k,1,j,i));
%                     yCompRotate_k=xComp(m,1,j,i).*sin((-1).*th_12_Z(k,1,j,i))+yComp(m,1,j,i).*cos(th_12_Z(k,1,j,i));
%                     
%                     
%                     
%                     
% %                     if(xCompRotate_j.*xCompRotate_k > 0.0)
% %                         if(yCompRotate_k > 0.0)
% %                             dihedralTheta(k,1,j,m,i)= (-1).*dihedralTheta(k,1,j,m,i);
% %                         end
% %                     end
% %                     
% %                     if(xCompRotate_j.*xCompRotate_k < 0.0)
% %                         if(yCompRotate_k < 0.0)
% %                             dihedralTheta(k,1,j,m,i)= (-1).*dihedralTheta(k,1,j,m,i);
% %                         end
% %                     end
%                     
%                     
%                     
%                     
%                     
%                     
%                     
%                     
%                     %                 if(i==501)
%                     %                     temp=0
%                     %                 end
%                     
%                     %                 if(theta_12_Z(j,1,2,i) < theta_12_Z(k,1,2,i))
%                     %                     dihedralTheta(j,1,2,k,i)= (-1).*dihedralTheta(j,1,2,k,i);
%                     %                 end
%                     
%                     
%                     
%                     
%                     
%                     
%                     
%                     %                 if(xComp(j,1,2,i) > 0.0 & xComp(k,1,2,i) > 0.0)%IF BOTH THE ATOMS j & k ARE IN 1st OR 4th QUADRANT IN NEWMAN PROJECTION
%                     %                     if(yComp(j,1,2,i) < yComp(k,1,2,i))
%                     %                         dihedralTheta(j,1,2,k,i)= -1.* dihedralTheta(j,1,2,k,i);
%                     %                     end
%                     %                 end
%                     %                 
%                     %                 if(xComp(j,1,2,i) < 0.0 & xComp(k,1,2,i) < 0.0)%IF BOTH THE ATOMS j & k ARE IN 2nd OR 3rd QUADRANT IN NEWMAN PROJECTION
%                     %                     if(yComp(j,1,2,i) > yComp(k,1,2,i))
%                     %                         dihedralTheta(j,1,2,k,i)= -1.* dihedralTheta(j,1,2,k,i);
%                     %                     end
%                     %                 end
%                     %                 
%                     %                 if(xComp(j,1,2,i) > 0.0 & xComp(k,1,2,i) < 0.0)%IF BOTH THE ATOMS j & k ARE IN 2nd OR 3rd QUADRANT IN NEWMAN PROJECTION
%                     % %                     if(yComp(j,1,2,i) > yComp(k,1,2,i))
%                     %                         dihedralTheta(j,1,2,k,i)= -1.* dihedralTheta(j,1,2,k,i);
%                     % %                     end
%                     %                 end
%                     %                 
%                     %                 if(yComp(j,1,2,i) > 0.0 & yComp(k,1,2,i) > 0.0)
%                     %                     if(xComp(j,1,2,i) > 0.0 & xComp(k,1,2,i) < 0.0)%IF THE ATOM j IS IN 1st QUADRANT & k IS IN 2nd QUADRANT IN NEWMAN PROJECTION
%                     %                         dihedralTheta(j,1,2,k,i)= -1.* dihedralTheta(j,1,2,k,i);
%                     %                     end
%                     %                 end
%                     %                 
%                     %                 if(yComp(j,1,2,i) < 0.0 & yComp(k,1,2,i) < 0.0)
%                     %                     if(xComp(j,1,2,i) < 0.0 & xComp(k,1,2,i) > 0.0)%IF THE ATOM j IS IN 3rd QUADRANT & k IS IN 4th QUADRANT IN NEWMAN PROJECTION
%                     %                         dihedralTheta(j,1,2,k,i)= -1.* dihedralTheta(j,1,2,k,i);
%                     %                     end
%                     %                 end
%                     
%                 end
%                 
%             end
%         end
%     end
% end 
% %ABOVE IS OLD CODE TO CALCULATE DIHEDRAL ANGLE

if( rem(fileNumMax,numConfig) ~= 0)
	fileNumMax = fileNumMax + 1;
end

configCount=1;
count=0;
lowerLimit=1;
upperLimit=numConfig;

for fileNum=1:fileNumMax
    file=strcat('si',num2str(fileNum));
    fileIn=strcat(file,'.com');
%     fileOut=strcat(file,'.log');
    
    fid=fopen(fileIn,'w');  % open a file to write
    % fprintf(fid,'%s','#!/bin/bash'); % unix file commands
    
    
    
    for i=lowerLimit:upperLimit
        
        fprintf(fid,'\n\n');
        fprintf(fid,'%s','g98 <<END > ','si_',num2str(configCount),'.log');  % specify a output file name for G98 with '.log' extension, used in G98 for output file
        fprintf(fid,'\n%s','%Mem=200MB');   %specify the memory
    
        fprintf(fid,'\n');
        fprintf(fid,'%s','%Chk=','si_',num2str(configCount));
        %fprintf(fid,'\n');
        %fprintf(fid,'%s','%RWF=/mnt/local1/si_',num2str(configCount));
		%         fprintf(fid,'\n');
		%         fprintf(fid,'%s','%SCR=/mnt/ramdisk/si_',num2str(configCount));
        fprintf(fid,'\n%s','%Nosave');

%         fprintf(fid,'\n%s','#p HF/STO-3G');
        fprintf(fid,'\n%s','#P B3LYP/6-31G** SCF=(Tight,MaxCycle=1024) Force');% specify the calculation level and the basis set
        fprintf(fid,'\n\n%s','silicon');% specify the title
        fprintf(fid,'\n\n%s','0 1');% specify the charge and the multiplicity
        
        fprintf(fid,'\n%s','Si1');% specify configurations in alternate z-matrix format for Gaussian 98
        
        fprintf(fid,'\n');
        fprintf(fid,'%s','Si2 Si1 ',r(1,2,i)); %this prints the number in the format 2.7000e00 rather than 2.70
        % r12=num2str(r(1,2,i)); % to print in the format 2.70, convert the number to string, which avoid printing like 2.7000e00 
        % fprintf(fid,'%s','Si2 Si1 ',r12); % print the string, DISADVANTAGE: IF r=90,THEN PRINTS 90 i.e.INTEGER AND NOT 90. or 90.0 WHICH IS REQUIRED IN Gaussian 98
        
        fprintf(fid,'\n');
        fprintf(fid,'%s','Si3 Si1 ',r(1,3,i),' Si2 ',theta(3,1,2,i));
        
        fprintf(fid,'\n');
        fprintf(fid,'%s','Si4 Si1 ',r(1,4,i),' Si2 ',theta(4,1,2,i),' Si3 ',dihedralTheta(4,1,2,3,i),' 0');%Z-matrix format using DIHEDRAL ANGLE
%         fprintf(fid,'%s','Si4 Si1 ',r(1,4,i),' Si2 ',theta(4,1,2,i),' Si3 ',theta(4,1,3,i),' 1'); % ALTERNATE z-MATRIX FORMAT
        
        fprintf(fid,'\n');
        fprintf(fid,'%s','Si5 Si1 ',r(1,5,i),' Si2 ',theta(5,1,2,i),' Si4 ',dihedralTheta(5,1,2,4,i),' 0'); %Z-matrix format using DIHEDRAL ANGLE
%         fprintf(fid,'%s','Si5 Si1 ',r(1,5,i),' Si2 ',theta(5,1,2,i),' Si4 ',theta(5,1,4,i),' 1');

       
  
%         fprintf(fid,'\n\n%s','--Link1--');
        %---------------------------------
        
%         fprintf(fid,'\n%s','%Mem=512MB');   %specify the memory
%         fprintf(fid,'\n');
%         fprintf(fid,'%s','%Chk=','si_',num2str(configCount));
%         fprintf(fid,'\n%s','#T UHF/6-31g** Guess=Read Geom=AllCheck SCF=(Tight,Maxcycle=128)');
%         
%         fprintf(fid,'\n\n%s','--Link1--');
        %---------------------------------

% %         fprintf(fid,'\n%s','%rwf=a,245mw,b,245mw,c,245mw,d,245mw,e,245mw,f,245mw,g,245mw,h,245mw');
%         fprintf(fid,'\n%s','%Mem=100MB');   %specify the memory
%         fprintf(fid,'\n');
%         fprintf(fid,'%s','%Chk=','si_',num2str(configCount));
%         fprintf(fid,'\n');
%         fprintf(fid,'%s','%RWF=/ramdisk/si_',num2str(configCount));
%         fprintf(fid,'\n%s','%Nosave');
%         fprintf(fid,'\n%s','#p MP4(SDQ)/6-31G** Guess=Read Geom=AllCheck SCF=(QC,TightLinEq) Force');% specify the calculation level and the basis set
        
%         if(i~=upperLimit)
            fprintf(fid,'\n\n%s','END');
            fprintf(fid,'\n');
            fprintf(fid,'%s','rm /local1/si_',num2str(configCount),'.rwf');
%             fprintf(fid,'\n');
%             fprintf(fid,'%s','rm /ramdisk/si_',num2str(configCount),'.scr');
%             fprintf(fid,'%s','rm si_',num2str(configCount),'.chk');
            fprintf(fid,'\n');
            fprintf(fid,'%s','echo "si_',num2str(configCount),'.log done"');
            fprintf(fid,'\n\n\n');
            
%         end%end if
        
%         if(i==upperLimit)
%             fprintf(fid,'\n\n%s','END');
%             fprintf(fid,'\n%s','echo "done"');
%         end%end if
        
        configCount=configCount+1;
        
    end%end for i
    
    
    
    % for i=1:total
    
    
    fclose(fid);
    
%     configCount=1;
    count=count+1;
    lowerLimit=(count.*numConfig)+1;
	
	if(lowerLimit ~= total-rem(total,numConfig)+1)
		upperLimit=lowerLimit+numConfig-1;
	else
		upperLimit=total;
	end
	
%     upperLimit=lowerLimit+numConfig-1;
    
    
    end % end for fileNum

