atoms=5;

training=200;
validate=0;
testing=0;


rThetaPhi=load('rThetaDihedralphi');

r=rThetaPhi.r;
theta=rThetaPhi.theta;
thetaDihedral=rThetaPhi.dihedralTheta;
% dihed=load('dihedralTh');
% thetaDihedral=dihed.dihedralTheta;


% [data,headertext]=xlsread('NNG98ConfigForceTRYWRONG.xls');

% config=data(:,1);

% r12=data(:,2);
% r13=data(:,3);
% r14=data(:,4);
% r15=data(:,5);

% r23=data(:,6);
% r24=data(:,7);
% r25=data(:,8);
% 
% r34=data(:,9);
% r35=data(:,10);
% 
% r45=data(:,11);

% th312=data(:,6);
% th412=data(:,7);
% th512=data(:,8);
% 
% phi4123=data(:,9);
% phi5124=data(:,10);





%ump4sdq=data(:,12);
% ump4sdq=data(:,22);
% 
% Fr12=data(:,23);
% Fr13=data(:,24);
% Fr14=data(:,25);
% Fr15=data(:,26);
% Fth312=data(:,27);
% Fth412=data(:,28);
% Fth512=data(:,29);
% Fphi4123=data(:,30);
% Fphi5124=data(:,31);


% f12=data(:,18);
% f13=data(:,19);
% f14=data(:,20);
% f15=data(:,21);
% 
% f12=f12';
% f13=f13';
% f14=f14';
% f15=f15';
% 
% f=[Fr12';Fr13';Fr14';Fr15';Fth312';Fth412';Fth512';Fphi4123';Fphi5124'];

temp=size(r);
total=temp(3);

% r=zeros(atoms,atoms,total);
% for i=1:total
%     r(1,2,i)=r12(i);
%     r(1,3,i)=r13(i);
%     r(1,4,i)=r14(i);
%     r(1,5,i)=r15(i);
%     
%     r(2,3,i)=r23(i);
%     r(2,4,i)=r24(i);
%     r(2,5,i)=r25(i);
%     
%     r(3,4,i)=r34(i);
%     r(3,5,i)=r35(i);
%     
%     r(4,5,i)=r45(i);
% end
% 
% for i=1:total
%     for j=1:atoms
%         for k=1:atoms
%             if(k~=j & k<j)
%                 r(j,k,i)=r(k,j,i);%Make the r matrix symmetric
%             end % end if
%         end % end for k
%     end % end for j
% end % end for i
% 
% 
% theta=zeros(atoms,atoms,atoms,total);
% for i=1:total
%     for j=1:atoms
%         for k=1:atoms
%             for m=1:atoms
%                 
%                 if (k~=j & m~=k & m~=j)
%                     %cosine(j,k,m,i)=(r(j,k,i).^2+r(k,m,i).^2-r(j,m,i).^2)./(2.*r(j,k,i).*r(k,m,i));
% %                     cosine(j,k,m,i)=(r(j,k,i).^2+r(j,m,i).^2-r(k,m,i).^2)./(2.*r(j,k,i).*r(j,m,i));
% %                     theta(j,k,m,i)=acos(cosine(j,k,m,i)).*180./pi;
%                     theta(j,k,m,i)=acos((r(j,k,i).^2+r(j,m,i).^2-r(k,m,i).^2)./(2.*r(j,k,i).*r(j,m,i))).*180./pi;
%                 end
%             end
%         end
%     end
% end

% ump4sdqConfig=ump4sdq';

%???????????????????????????
% totDEL=0;
% for i=1:total
% 	if(ump4sdqConfig(i) == 0)
% 		totDEL=totDEL+1;
% 		delet(totDEL)=i;
% 	end
% end
%???????????????????????????

% total=total-totDEL;
% 
% 
% for i=1:total
% 	for m=1:totDel
% 		if(i == delet(m))
			
		
		

rr1=zeros((atoms-1),total);

for i=1:total
    m=1;
%     for j=1
        for k=2:atoms
            rr1(m,i)=r(1,k,i);
            
            m=m+1;
        end
%     end
i
end
    
order_sort=zeros(4,total);
for i=1:total
    order_sort(1,i)=2;
    order_sort(2,i)=3;
    order_sort(3,i)=4;
    order_sort(4,i)=5;
end



% rThConfig=zeros((3*atoms-6),total);
% 
% for i=1:total
%     m=1;
%     for k=1:(atoms-1)
%         rThConfig(m,i)=rr1(k,i);
%         m=m+1;
%     end
% 	i
% end
% 
% for i=1:total
%     rThConfig(5,i)=theta(2,1,3,i);
%     rThConfig(6,i)=theta(2,1,4,i);
%     rThConfig(7,i)=theta(2,1,5,i);
%     rThConfig(8,i)=thetaDihedral(4,1,2,3,i);
%     rThConfig(9,i)=thetaDihedral(5,1,2,4,i);
% end
% 
% rTh_ump4sdqConfig_UNSORT=[rThConfig ; ump4sdqConfig];
% 
% for i=1:total
% 	if(rTh_ump4sdqConfig_UNSORT(8,i) < 0)
% 		rTh_ump4sdqConfig_UNSORT(8,i) = -1.*rTh_ump4sdqConfig_UNSORT(8,i);
% 		rTh_ump4sdqConfig_UNSORT(9,i) = -1.*rTh_ump4sdqConfig_UNSORT(9,i);
% 	end
% end


fConfig=zeros((atoms-1),total);






rr1_sort=rr1;

% for i=1:total
%     for j=1:(atoms-1)-1
%         for k=1:(atoms-1)-j
%             if(rr1_sort(k,i) > rr1_sort(k+1,i))
%                 temp=rr1_sort(k,i);
%                 rr1_sort(k,i)=rr1_sort(k+1,i);
%                 rr1_sort(k+1,i)=temp;
%                 
%                 temp_order=order_sort(k,i);
%                 order_sort(k,i)=order_sort(k+1,i);
%                 order_sort(k+1,i)=temp_order;
%                 
%             end
%         end
%     end
% end

rThConfig_sort=zeros((3*atoms-6),total);

for i=1:total
    m=1;
%     for j=1:(atoms-1)
        for k=1:(atoms-1)
            rThConfig_sort(m,i)=rr1_sort(k,i);
            m=m+1;
        end
%     end
end
  
for i=1:total
    rThConfig_sort(5,i)=theta(order_sort(2,i),1,order_sort(1,i),i);
    rThConfig_sort(6,i)=theta(order_sort(3,i),1,order_sort(1,i),i);
    rThConfig_sort(7,i)=theta(order_sort(4,i),1,order_sort(1,i),i);
    rThConfig_sort(8,i)=thetaDihedral(order_sort(3,i),1,order_sort(1,i),order_sort(2,i),i);
    rThConfig_sort(9,i)=thetaDihedral(order_sort(4,i),1,order_sort(1,i),order_sort(3,i),i);
end
    
for i=1:total
	if(rThConfig_sort(8,i) < 0)
		rThConfig_sort(8,i) = -1.*rThConfig_sort(8,i);
		rThConfig_sort(9,i) = -1.*rThConfig_sort(9,i);
		f(8,i) = -1.*f(8,i);
		f(9,i) = -1.*f(9,i);
	end
end




% rThConfig_sort_dihedralSign=zeros((3*atoms-6),2*total);
% ump4sdqConfig_sort_dihedralSign=zeros(1,2*total);
% 
% for i=1:total
% 
%         rThConfig_sort_dihedralSign(1,2*i-1)=rThConfig_sort(1,i);
%         rThConfig_sort_dihedralSign(2,2*i-1)=rThConfig_sort(2,i);
%         rThConfig_sort_dihedralSign(3,2*i-1)=rThConfig_sort(3,i);
%         rThConfig_sort_dihedralSign(4,2*i-1)=rThConfig_sort(4,i);
%         rThConfig_sort_dihedralSign(5,2*i-1)=rThConfig_sort(5,i);
%         rThConfig_sort_dihedralSign(6,2*i-1)=rThConfig_sort(6,i);
%         rThConfig_sort_dihedralSign(7,2*i-1)=rThConfig_sort(7,i);
%         rThConfig_sort_dihedralSign(8,2*i-1)=rThConfig_sort(8,i);
%         rThConfig_sort_dihedralSign(9,2*i-1)=rThConfig_sort(9,i);
% 
%         ump4sdqConfig_sort_dihedralSign(2*i-1)=ump4sdqConfig(i);
%         
%         rThConfig_sort_dihedralSign(1,2*i)=rThConfig_sort(1,i);
%         rThConfig_sort_dihedralSign(2,2*i)=rThConfig_sort(2,i);
%         rThConfig_sort_dihedralSign(3,2*i)=rThConfig_sort(3,i);
%         rThConfig_sort_dihedralSign(4,2*i)=rThConfig_sort(4,i);
%         rThConfig_sort_dihedralSign(5,2*i)=rThConfig_sort(5,i);
%         rThConfig_sort_dihedralSign(6,2*i)=rThConfig_sort(6,i);
%         rThConfig_sort_dihedralSign(7,2*i)=rThConfig_sort(7,i);
%         rThConfig_sort_dihedralSign(8,2*i)= -1.* rThConfig_sort(8,i);
%         rThConfig_sort_dihedralSign(9,2*i)= -1.* rThConfig_sort(9,i);
%         
%         ump4sdqConfig_sort_dihedralSign(2*i)=ump4sdqConfig(i);
% 
% end
        
        


% fConfig_sort=zeros((atoms-1),total);
% 
% for i=1:total
%     for k=1:(atoms-1)
%         fConfig_sort(k,i)=f(order_sort(k,i)-1,i);
%     end
% end

    
% f1=f(1,:);


% r_order=zeros(atoms,atoms,factorial(atoms-1),total);
% 
% f_order=zeros((atoms-1),factorial(atoms-1),total);
% 
% ump4sdq_order=zeros(1,factorial(atoms-1),total);
% 
% for i=1:total
%     
%     ordering=1;
%     
%     %     for ordering=1:factorial(atoms-1)
%     while(ordering<=factorial(atoms-1))
%         for order2=2:atoms
%             for order3=2:atoms
%                 if(order3 ~= order2)
%                     for order4=2:atoms
%                         if(order4 ~= order3 & order4 ~= order2)
%                             for order5=2:atoms
%                                 if(order5 ~= order4 & order5 ~= order3 & order5 ~= order2)
%                                     
%                                     r_order(1,2,ordering,i)=r(1,order2,i);
%                                     r_order(1,3,ordering,i)=r(1,order3,i);
%                                     r_order(1,4,ordering,i)=r(1,order4,i);
%                                     r_order(1,5,ordering,i)=r(1,order5,i);
%                                     
%                                     r_order(2,3,ordering,i)=r(order2,order3,i);
%                                     r_order(2,4,ordering,i)=r(order2,order4,i);
%                                     r_order(2,5,ordering,i)=r(order2,order5,i);
%                                     
%                                     r_order(3,4,ordering,i)=r(order3,order4,i);
%                                     r_order(3,5,ordering,i)=r(order3,order5,i);
%                                     
%                                     r_order(4,5,ordering,i)=r(order4,order5,i);
%                                     
%                                     
%                                     f_order(1,ordering,i)=f((order2-1),i);
%                                     f_order(2,ordering,i)=f((order3-1),i);
%                                     f_order(3,ordering,i)=f((order4-1),i);
%                                     f_order(4,ordering,i)=f((order5-1),i);
%                                     
%                                     ump4sdq_order(1,ordering,i)=ump4sdq(i); 
% 	                             
%                                     ordering=ordering+1;
%                                     
%                                 end% end if
%                             end% end for order5
%                         end% end if 
%                     end% end for order4
%                 end% end if(oredr3 ~= order2) 
%             end% end for order3
%         end% end for order2
%     end% end while ordering   
% end% end for i




% rTh=zeros(((atoms*(atoms-1)./2)-1),factorial(atoms-1),total);
% 
% for i=1:total
%     
%     ordering=1;
%     
%     %     for ordering=1:factorial(atoms-1)
%     while(ordering<=factorial(atoms-1))
%         for order2=2:atoms
%             for order3=2:atoms
%                 if(order3 ~= order2)
%                     for order4=2:atoms
%                         if(order4 ~= order3 & order4 ~= order2)
%                             for order5=2:atoms
%                                 if(order5 ~= order4 & order5 ~= order3 & order5 ~= order2)
%                                     
%                                     rTh(1,ordering,i)=r(1,order2,i);
%                                     rTh(2,ordering,i)=r(1,order3,i);
%                                     rTh(3,ordering,i)=r(1,order4,i);
%                                     rTh(4,ordering,i)=r(1,order5,i);
%                                     
%                                     rTh(5,ordering,i)=theta(1,order2,order3,i);
%                                     rTh(6,ordering,i)=theta(1,order2,order4,i);
%                                     rTh(7,ordering,i)=theta(1,order2,order5,i);
%                                     rTh(8,ordering,i)=theta(1,order3,order4,i);
%                                     rTh(9,ordering,i)=theta(1,order4,order5,i);
%                                     
%                                     ordering=ordering+1;
%                                     
%                                 end% end if
%                             end% end for order5
%                         end% end if 
%                     end% end for order4
%                 end% end if(oredr3 ~= order2) 
%             end% end for order3
%         end% end for order2
%     end% end while ordering   
% end% end for i




% rThConfig=zeros(((atoms*(atoms-1)./2)-1),total*factorial(atoms-1));
% 
% m=1;
% for i=1:total
%     for ordering=1:factorial(atoms-1)
%         for j=1:((atoms*(atoms-1)./2)-1)
%             rThConfig(j,m)=rTh(j,ordering,i);
%         end
%         m=m+1;
%     end
% end


% rConfig=zeros(atoms,atoms,total*factorial(atoms-1));
% 
% m=1;
% for i=1:total
%     
%     for ordering=1:factorial(atoms-1)
%         for j=1:atoms
%             for k=1:atoms
%                 rConfig(j,k,m)=r_order(j,k,ordering,i);
%                 
%             end
%         end
% 
%         m=m+1
%     end
%     
% end


% fConfig=zeros(atoms-1,total*factorial(atoms-1));
% 
% m=1;
% for i=1:total
%     for ordering=1:factorial(atoms-1)
%         for j=1:atoms-1
%             fConfig(j,m)=f_order(j,ordering,i);
%         end% end for j
%         m=m+1
%     end%end for ordering
% end%end for i
% 
% f1=zeros(1,total*factorial(atoms-1));
% 
% for i=1:total*factorial(atoms-1)
%     f1(1,i)=fConfig(1,i);
% end
% 
% 
% 
% ump4sdqConfig=zeros(1,total*factorial(atoms-1));
% 
% m=1;
% for i=1:total
%     for ordering=1:factorial(atoms-1)
%         ump4sdqConfig(1,m)=ump4sdq_order(1,ordering,i);
%         m=m+1;
%     end
% end




%  r1=zeros(atoms*(atoms-1)/2-1,total*factorial(atoms-1));
% %r1=zeros(atoms*(atoms-1)/2,total*factorial(atoms-1));
% for i=1:total*factorial(atoms-1)
%     m=1;
%     for j=1:atoms-2
%     %for j=1:atoms
%         for k=1:atoms
%             %	for m=1:(atoms*atoms-atoms)/2
%             
%             if (k~=j & k>j)
%                 r1(m,i)=(rConfig(j,k,i));%r1(m,i) STORES ALL n(n-1)/2 BOND DISTANCES
%                 m=m+1;
%             end % end if
%             
%             i            
%         end
%     end
% end




                                    

 %[r_in,minr,maxr,f_in,minf,maxf]=premnmx(r1,ump4sdqConfig);

%  [r_in,minr,maxr,f_in,minf,maxf]=premnmx(rThConfig_sort,fConfig_sort);
% % %r_in=r1;
% % % v_in=v1;
% % 
% % dv1dr_train=dv1dr_in(:,1:training);
% % dv1dr_test=dv1dr_in(:,training+validate+testing+1:total);
% % 
% r_train=r_in(:,1:training);
% val.P=r_in(:,training+1:training+validate);
% test.P=r_in(:,training+validate+1:training+validate+testing);
% r_test=r_in(:,training+validate+testing+1:total);
% % % 
% f_train=f_in(:,1:training);
% val.T=f_in(:,training+1:training+validate);
% test.T=f_in(:,training+validate+1:training+validate+testing);
% f_test=f_in(:,training+validate+testing+1:total);



% clear r12; clear r13; clear r14; clear r15; clear r23; clear r24; clear r25; clear r34; clear r35; clear r45;
% clear f12; clear f13; clear f14; clear f15;
%clear r; clear r_order; clear rConfig;
%clear f, clear f_order;

% f12_order=zeros(total*factorial(atoms-1));
% f13_order=zeros(total*factorial(atoms-1));
% f14_order=zeros(total*factorial(atoms-1));
% f15_order=zeros(total*factorial(atoms-1));

% m=1;
% for i=1:total*factorial(atoms-1)
% 
%     for m=(i-1)*24+1:i*24 
%         f12_order(m)=f12(i);
%     end
%     
%     for m=(i-1)*24+1:i*24 
%         f13_order(m)=f13(i);
%     end
%     
%     for m=(i-1)*24+1:i*24 
%         f14_order(m)=f14(i);
%     end
%     
%     for m=(i-1)*24+1:i*24 
%         f15_order(m)=f15(i);
%     end
%     
%     i
% end
% 
% fConfig=[f12_order;f13_order;f14_order;f15_order];
        



% rTh_ump4sdqConfig=[rThConfig_sort_dihedralSign ; ump4sdqConfig_sort_dihedralSign];
rTh_ump4sdqConfig=[rThConfig_sort ; ump4sdqConfig];

%???????????????????????????
rTh_ump4sdqConfig=[rTh_ump4sdqConfig(:,1:331) rTh_ump4sdqConfig(:,333:2095) rTh_ump4sdqConfig(:,2097:9216) rTh_ump4sdqConfig(:,9218:total)];
%???????????????????????????

fid_config=fopen('NNG98ConfigSORTPermute.xls','w');

fprintf(fid_config,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s','r12','r13','r14','r15','th3-1-2','th4-1-2','th5-1-2','phi4-1-2-3','phi5-1-2-4','ump4(sdq)');
fprintf(fid_config,'\n');

for i=1:(total)
    for j=1:(3*atoms-6)+1
        fprintf(fid_config,'%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f',rTh_ump4sdqConfig(j,i));
    end
    fprintf(fid_config,'\n');
end


fclose(fid_config);
