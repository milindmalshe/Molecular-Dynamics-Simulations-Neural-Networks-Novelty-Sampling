
function [coord,vel,accelrn,tersoff_PE,KE]= Gear_predicorrectSi5(count,coord,total,numMov,numPeriph,numBound,movAtom,periphAtom,boundAtom,vel,accelrn,mass);

global  countTotal;
global B_accelrn C_accelrn;
%  if(count <= 250) 
 	delT = 0.1;
%  else
%  	delT = -0.05;
%  end

KE = 0;
% accelrn = zeros(total,3);%????????????????????????
% [tersoff_PE,force] = NNG98_2(coord,total,numMov,numPeriph,numBound,movAtom,periphAtom,boundAtom);
% 
% accelrn = -1.*force./mass;


if(count == 1)
	B_accelrn = zeros(total,3);
	C_accelrn = zeros(total,3);
end

c1=delT;
c2=c1*delT/2;
c3=c2*delT/3;
c4=c3*delT/4;

GEAR0=19.0 / 120.0; GEAR1=3.0 / 4.0; GEAR3=1.0 / 2.0; GEAR4 = 1.0 / 12.0;
CR = GEAR0 * c2;
CV = GEAR1 * c2 / c1;
CB = GEAR3 * c2 / c3;
CC = GEAR4 * c2 / c4;

coord = coord + c1*vel + c2*accelrn + c3*B_accelrn + c4*C_accelrn;
vel = vel + c1*accelrn + c2*B_accelrn + c3*C_accelrn;
accelrn = accelrn + c1*B_accelrn + c2*C_accelrn;
B_accelrn = B_accelrn + c1*C_accelrn;

           
% coord= coord + vel.*delT + 1/2.*accelrn.*delT^2;

% vel= vel + 1./2. *accelrn.*(delT);

% [coord,a,f,pe] = forceLJ(coord,a,f,pe);

[tersoff_PE,force] = tersoffSi3(coord,total,numMov,numPeriph,numBound,movAtom,periphAtom,boundAtom);

accelrn_NEW = -1.*force./mass;
corrector = accelrn_NEW - accelrn;

coord = coord + CR .* corrector;
vel = vel + CV .* corrector;
accelrn = accelrn_NEW;
B_accelrn = B_accelrn + CB * corrector;
C_accelrn = C_accelrn + CC * corrector;

% vel= vel + 1./2. *accelrn.*(delT);

for i=1:total
	velTot(i)= sqrt(sum(vel(i,:).^2));
end

for i=1:total
	KE= KE + 1/2*mass.*(velTot(i).^2);
end

% KE= 1./2.* mass.* sum(vel.^2);
KE;

