

coordStoreTersoff=load('coordStore_tersoffOLD.mat');
coordStore=coordStoreTersoff.coordStore;

fidx=fopen('x5_coordStore_tersoffOLD.xls','w');
fidy=fopen('y5_coordStore_tersoffOLD.xls','w');
fidz=fopen('z5_coordStore_tersoffOLD.xls','w');


temp = size(coordStore);
total = temp(3);

for i= 1:total
% 	for j=2:5
		if(rem(i,1)==0)

		x2 = coordStore(2,1,i)-coordStore(1,1,i);
		y2 = coordStore(2,2,i)-coordStore(1,2,i);
		z2 = coordStore(2,3,i)-coordStore(1,3,i);
		
		x3 = coordStore(3,1,i)-coordStore(1,1,i);
		y3 = coordStore(3,2,i)-coordStore(1,2,i);
		z3 = coordStore(3,3,i)-coordStore(1,3,i);
		
		x4 = coordStore(4,1,i)-coordStore(1,1,i);
		y4 = coordStore(4,2,i)-coordStore(1,2,i);
		z4 = coordStore(4,3,i)-coordStore(1,3,i);
		
		x5 = coordStore(5,1,i)-coordStore(1,1,i);
		y5 = coordStore(5,2,i)-coordStore(1,2,i);
		z5 = coordStore(5,3,i)-coordStore(1,3,i);
		
		
		fprintf(fidx,'%f\t%f\t%f\t%f', x2 , x3 , x4 , x5);
		% 	fprintf(fidx,'\n'); fprintf(fidy,'\n'); fprintf(fidz,'\n'); 
		fprintf(fidy,'%f\t%f\t%f\t%f', y2 , y3 , y4 , y5);
		% 	fprintf(fidx,'\n'); fprintf(fidy,'\n'); fprintf(fidz,'\n'); 
		fprintf(fidz,'%f\t%f\t%f\t%f', z2 , z3 , z4 , z5);
		fprintf(fidx,'\n'); 
		fprintf(fidy,'\n');
		fprintf(fidz,'\n'); 
		
		i
	end
% 	end
end

fclose(fidx);
fclose(fidy);
fclose(fidz);