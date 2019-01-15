% function [coord,v,a,f,pe,ke]= velverlet(coord,v,a,f,pe,ke);

function [coord,vel,accelrnNEW,tersoff_PE,KE]= velverlet_compPhys(count,coord,total,numMov,numPeriph,numBound,numSurface,movAtom,periphAtom,boundAtom,surfaceAtom,vel,accelrn,mass);

global delT countTotal;
% if(count <= 250) 
	delT = 0.1;
% else
% 	delT = -0.05;
% end

KE = 0;
% accelrn = zeros(total,3);%????????????????????????

coord= coord + vel.*delT + 1/2.*accelrn.*delT^2;

% vel= vel + 1./2. *accelrn.*(delT);

% [coord,a,f,pe] = forceLJ(coord,a,f,pe);


[tersoff_PE,force] = NNG98_2_coordStore_PBC_4ForceComponentOnALLatomsInList(coord,total,numMov,numPeriph,numBound,numSurface,movAtom,periphAtom,boundAtom,surfaceAtom);
% [tersoff_PE,force] = tersoffSi3_PeriodicBoundaryCond(coord,total,numMov,numPeriph,numBound,numSurface,movAtom,periphAtom,boundAtom,surfaceAtom);
% [tersoff_PE,force] = NNG98_2_coordStore(coord,total,numMov,numPeriph,numBound,movAtom,periphAtom,boundAtom);

accelrnNEW = -1.*force./mass;

vel= vel + 1./2. *(accelrn + accelrnNEW).*(delT);

for i=1:total
	velTot(i)= sqrt(sum(vel(i,:).^2));
end

for i=1:total
	KE= KE + 1/2*mass.*(velTot(i).^2);
end

% KE= 1./2.* mass.* sum(vel.^2);
KE;




