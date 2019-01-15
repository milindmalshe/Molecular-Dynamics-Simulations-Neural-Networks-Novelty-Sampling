function [list] = bondList(coord,numMov,total,rCutoff,movAtom)

% global x y z;
x=coord(:,1);
y=coord(:,2);
z=coord(:,3);

for i=1:numMov
	iMov=movAtom(i);
	numBond=0;
	for j=1:total
		if(j ~= iMov)
			r= sqrt( (x(iMov)-x(j)).^2 + (y(iMov)-y(j)).^2 + (z(iMov)-z(j)).^2 );
			r;
			if(r <= (rCutoff))
				numBond=numBond+1;
				list(iMov,1)=numBond;
				list(iMov,numBond+1)=j;
			end
		end
	end
end