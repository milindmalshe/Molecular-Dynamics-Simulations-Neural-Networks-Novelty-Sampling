
function [Vij,DEDmov] = tersoffSi3_OLD(coord,total,numMov,numPeriph,numBound,movAtom,periphAtom,boundAtom);

% global x y z;
global A  B lamda1 lamda2 lamda3 alpha beta eta c d h R D rCutoff;
global DrijDxi DrijDyi DrijDzi DrikDxi DrikDyi DrikDzi DrjkDxj DrjkDyj DrjkDzj; 

global countStore;% ****THIS MODIFICATION IS MADE WHILE STORING THE CONFIGURATIONS 
global coordStore;% ****THIS MODIFICATION IS MADE WHILE STORING THE CONFIGURATIONS 
global countFFT;
global rijFFT;

A=3.2647e3;
B=9.5373e1;
lamda1=3.2394;
lamda2=1.3258;
lamda3=1.3258;
alpha=0;
beta=3.3675e-1;
eta=2.2956e1;
c=4.8381;
d=2.0417;
h=0;
R=3.0;
D=0.2; 

rCutoff=3.0;

Vij=0;


% for i=1:total
% 	DEDx(i)=0;
% 	DEDy(i)=0;
% 	DEDz(i)=0;
% 	
% 	Dzetak= zeros(total,3);
% end

DED= zeros(total,3);
Dzetak= zeros(total,3);

Dbij=zeros(2,3);
Dbk= zeros(total,3);

x=coord(:,1);
y=coord(:,2);
z=coord(:,3);

list=bondList(coord,numMov,total,rCutoff,movAtom); 

% countStore = 0; % ****THIS MODIFICATION IS MADE WHILE STORING THE CONFIGURATIONS 

for i=1:numMov
	iMov=movAtom(i);
% 	if(list(iMov,1) ~=1 & list(iMov,1) ~=2 & list(iMov,1) ~=4)
% 		iMov
% 	end
	for j=2:list(iMov,1)+1
		jBond=list(iMov,j);
		if(jBond == iMov)
			continue;
		end
		rij= sqrt( (x(iMov)-x(jBond)).^2 + (y(iMov)-y(jBond)).^2 + (z(iMov)-z(jBond)).^2 );
		
		DrijDxi= (x(iMov)-x(jBond))./rij;
		DrijDyi= (y(iMov)-y(jBond))./rij;
		DrijDzi= (z(iMov)-z(jBond))./rij;
		
		[fR,Dfrij] = fr(rij);
		[fA,Dfaij] = fa(rij);
		% 			[fC,Dfcij] = fc(rij,iMov,jBond);
		
		
		zetaij=0;
		Dzetaij=zeros(2,3);
		% 			Dzetak= zeros(total,3);
		
		%START 3-body parameters calculation
		for k=2:list(iMov,1)+1
			k3Body=list(iMov,k);
			if((k3Body == iMov) | (k3Body == jBond))
				continue
			end
			rik= sqrt( (x(iMov)-x(k3Body)).^2 + (y(iMov)-y(k3Body)).^2 + (z(iMov)-z(k3Body)).^2 );
			rjk= sqrt( (x(jBond)-x(k3Body)).^2 + (y(jBond)-y(k3Body)).^2 + (z(jBond)-z(k3Body)).^2 );
			
			
			DrikDxi= (x(iMov)-x(k3Body))./rik;
			DrikDyi= (y(iMov)-y(k3Body))./rik;
			DrikDzi= (z(iMov)-z(k3Body))./rik;
			
			DrjkDxj= (x(jBond)-x(k3Body))./rjk;%????????????????????????????
			DrjkDyj= (y(jBond)-y(k3Body))./rjk;
			DrjkDzj= (z(jBond)-z(k3Body))./rjk;
			
			[zetaij,Dzetaij,Dzetak] = zeta(k3Body,rij,rik,rjk,zetaij,Dzetaij,Dzetak);		
			
% 			[fC,Dfcik] = fc(rik,iMov,k3Body);
% 			
% 			[g_ijk,Dg] = g(rij,rik,rjk);
% 			
% 			zeta=zeta + fC*g_ijk;
% 			
% 			Dzetaij(1,1) = Dzetaij(1,1) + (fC*Dg(1,1) + Dfcik(1,1)*g_ijk);%derivative w.r.r. xi
% 			Dzetaij(1,2) = Dzetaij(1,2) + (fC*Dg(1,2) + Dfcik(1,2)*g_ijk);
% 			Dzetaij(1,3) = Dzetaij(1,3) + (fC*Dg(1,3) + Dfcik(1,3)*g_ijk);
% 			
% 			Dzetaij(2,1) = Dzetaij(2,1) + fC*Dg(2,1);%derivative w.r.r. xj
% 			Dzetaij(2,2) = Dzetaij(2,2) + fC*Dg(2,2);
% 			Dzetaij(2,3) = Dzetaij(2,3) + fC*Dg(2,3);
% 			
% 			Dzetak(k3Body,1)= fC*Dg(3,1) + g_ijk*Dfcik(2,1);%derivative w.r.r. xk
% 			Dzetak(k3Body,2)= fC*Dg(3,2) + g_ijk*Dfcik(2,2);
% 			Dzetak(k3Body,3)= fC*Dg(3,3) + g_ijk*Dfcik(2,3);
% 			
% 			% 				end
		end
		
		

		[bij,Dbij,Dbk] = b(iMov,jBond,zetaij,Dzetaij,Dzetak,list,Dbij,Dbk);
		
% 		bij=(1+ (beta^eta * zetaij^eta))^(-1/(2*eta));
% 		
% 		if(zetaij == 0)
% 			zetaij = 1.0e-10; % to avoid infinity if zeta is 0 then in the calculation of derivative of the bij index of beta becomes -ve, therefore it becomes infinity
% 		end
% 		
% 		% 			if(zeta ~= 0)
% 		Dbij(1,1) = -1/2*((1+(beta^eta)*(zetaij^eta))^(-1/2*(1+2*eta)/eta))*(beta^eta)*zetaij^(eta-1) .* Dzetaij(1,1);
% 		Dbij(1,2) = -1/2*((1+(beta^eta)*(zetaij^eta))^(-1/2*(1+2*eta)/eta))*(beta^eta)*zetaij^(eta-1) .* Dzetaij(1,2);
% 		Dbij(1,3) = -1/2*((1+(beta^eta)*(zetaij^eta))^(-1/2*(1+2*eta)/eta))*(beta^eta)*zetaij^(eta-1) .* Dzetaij(1,3);
% 		
% 		Dbij(2,1) = -1/2*((1+(beta^eta)*(zetaij^eta))^(-1/2*(1+2*eta)/eta))*(beta^eta)*zetaij^(eta-1) .* Dzetaij(2,1);
% 		Dbij(2,2) = -1/2*((1+(beta^eta)*(zetaij^eta))^(-1/2*(1+2*eta)/eta))*(beta^eta)*zetaij^(eta-1) .* Dzetaij(2,2);
% 		Dbij(2,3) = -1/2*((1+(beta^eta)*(zetaij^eta))^(-1/2*(1+2*eta)/eta))*(beta^eta)*zetaij^(eta-1) .* Dzetaij(2,3);
% 		
% 		for k=2:list(iMov,1)+1
% 			k3Body=list(iMov,k);
% 			if((k3Body == iMov) | (k3Body == jBond))
% 				continue;
% 			end
% 			Dbk(k3Body,1) = -1/2*((1+(beta^eta)*(zetaij^eta))^(-1/2*(1+2*eta)/eta))*(beta^eta)*zetaij^(eta-1) .* Dzetak(k3Body,1);
% 			Dbk(k3Body,2) = -1/2*((1+(beta^eta)*(zetaij^eta))^(-1/2*(1+2*eta)/eta))*(beta^eta)*zetaij^(eta-1) .* Dzetak(k3Body,2);
% 			Dbk(k3Body,3) = -1/2*((1+(beta^eta)*(zetaij^eta))^(-1/2*(1+2*eta)/eta))*(beta^eta)*zetaij^(eta-1) .* Dzetak(k3Body,3);
% 			% 				end
% 		end
		
		% 			else
		% 				Dbij(1,1)=0;
		% 				Dbij(1,2)=0;
		% 				Dbij(1,3)=0;
		% 				
		% 				Dbij(2,1)=0;
		% 				Dbij(2,2)=0;
		% 				Dbij(2,3)=0;
		% 				
		% 				Dbk(k3Body,1)=0;
		% 				Dbk(k3Body,2)=0;
		% 				Dbk(k3Body,3)=0;
		% 			end
		
		
		
		
		[fC,Dfcij] = fc(rij);
		
		Vij = Vij + fC*(fR + bij * fA); 
		
		DED(iMov,1) = DED(iMov,1)+ Dfcij(1,1).*(fR + bij*fA) + (fC.*(Dfrij(1,1)+(bij.*Dfaij(1,1))+(fA.*Dbij(1,1))));
		DED(iMov,2) = DED(iMov,2)+ Dfcij(1,2).*(fR + bij*fA) + (fC.*(Dfrij(1,2)+(bij.*Dfaij(1,2))+(fA.*Dbij(1,2))));
		DED(iMov,3) = DED(iMov,3)+ Dfcij(1,3).*(fR + bij*fA) + (fC.*(Dfrij(1,3)+(bij.*Dfaij(1,3))+(fA.*Dbij(1,3))));
		
		DED(jBond,1) = DED(jBond,1)+ Dfcij(2,1).*(fR + bij*fA) + (fC.*(Dfrij(2,1)+(bij.*Dfaij(2,1))+(fA.*Dbij(2,1))));
		DED(jBond,2) = DED(jBond,2)+ Dfcij(2,2).*(fR + bij*fA) + (fC.*(Dfrij(2,2)+(bij.*Dfaij(2,2))+(fA.*Dbij(2,2))));
		DED(jBond,3) = DED(jBond,3)+ Dfcij(2,3).*(fR + bij*fA) + (fC.*(Dfrij(2,3)+(bij.*Dfaij(2,3))+(fA.*Dbij(2,3))));
		
		for k=2:list(iMov,1)+1
			k3Body=list(iMov,k);
			if((k3Body == iMov) | (k3Body == jBond))
				continue;
			end
			
			DED(k3Body,1) = DED(k3Body,1)+ fC.*Dbk(k3Body,1).*fA;
			DED(k3Body,2) = DED(k3Body,2)+ fC.*Dbk(k3Body,2).*fA;
			DED(k3Body,3) = DED(k3Body,3)+ fC.*Dbk(k3Body,3).*fA;
			% 				end
		end
		
		
		% 		end
	end   %end for jBond
	
	% % ****FOLLOWING THIS MODIFICATION IS MADE WHILE STORING THE CONFIGURATIONS 
% 	if(iMov==14 | iMov==44 | iMov==47 | iMov==56 | iMov==59 | iMov==84 | iMov==85 | iMov==90 | iMov==91 | iMov==110 | iMov==111 | iMov==112 | iMov==113 | iMov==137 | iMov==164 | iMov==191 | iMov==218)
	if(iMov==1 | iMov==71 | iMov==73 | iMov==82 | iMov==92)
		if(list(iMov,1) == 4)
			countStore = countStore+1;
			for j=2:list(iMov,1)+1
				jBond=list(iMov,j);
				if(jBond == iMov)
					continue;
				end
				rij= sqrt( (x(iMov)-x(jBond)).^2 + (y(iMov)-y(jBond)).^2 + (z(iMov)-z(jBond)).^2 );		
				% 				for k=2:list(iMov,1)+1
				% 					k3Body=list(iMov,k);
				% 					if((k3Body == iMov) | (k3Body == jBond))
				% 						continue;
				% 					end
				
				coordStore(1,1,countStore) = x(iMov); coordStore(1,2,countStore) = y(iMov); coordStore(1,3,countStore) = z(iMov);
				coordStore(j,1,countStore) = x(jBond); coordStore(j,2,countStore) = y(jBond); coordStore(j,3,countStore) = z(jBond);
				% 					coordStore(3,1,countStore) = x(jBond); coordStore(3,2,countStore) = y(jBond); coordStore(3,3,countStore) = z(jBond);
				% 					coordStore(4,1,countStore) = x(jBond); coordStore(4,2,countStore) = y(jBond); coordStore(4,3,countStore) = z(jBond);
				% 					coordStore(5,1,countStore) = x(jBond); coordStore(5,2,countStore) = y(jBond); coordStore(5,3,countStore) = z(jBond);
				% 					coordStore(6,1,countStore) = Vij;
                % 				end %end for k
            end% end for j
        end
    end % end if(iMov== | iMov==....)
    % % ****ABOVE THIS MODIFICATION IS MADE WHILE STORING THE CONFIGURATIONS 
    
    if(iMov ==14)
        countFFT = countFFT+1;
        for j=2:list(iMov,1)+1
            jBond=list(iMov,j);
            if(jBond == iMov)
                continue;
            end
            %             rij= sqrt( (x(iMov)-x(jBond)).^2 + (y(iMov)-y(jBond)).^2 + (z(iMov)-z(jBond)).^2 );	
            rijFFT(countFFT,j-1)= sqrt( (x(iMov)-x(jBond)).^2 + (y(iMov)-y(jBond)).^2 + (z(iMov)-z(jBond)).^2 );                           
        end
        
        j=2;
        jBond=list(iMov,j);
        for k=3:list(iMov,1)+1
            k3Body=list(iMov,k);
            if((k3Body == iMov) | (k3Body == jBond))
                continue;
            end
            rjk = sqrt( (x(jBond)-x(k3Body)).^2 + (y(jBond)-y(k3Body)).^2 + (z(jBond)-z(k3Body)).^2 );
            rijFFT(countFFT,(list(iMov,1)+(k-2))) = acos((rijFFT(j-1).^2 + rijFFT(k-1).^2 - rjk.^2) ./ (2.* rijFFT(j-1) .* rijFFT(k-1))).*180./pi;% consine rule 
        end
    end
end   %end for iMov

Vij = Vij/2;

for i=1:total
	DED(i,1) = DED(i,1)./2;
	DED(i,2) = DED(i,2)./2;
	DED(i,3) = DED(i,3)./2;
end

DEDmov = zeros(total,3);
for i=1:numMov
	iMov=movAtom(i);
	DEDmov(iMov,1) = DED(iMov,1);
	DEDmov(iMov,2) = DED(iMov,2);
	DEDmov(iMov,3) = DED(iMov,3);
end

% DE=[DEDx' DEDy' DEDz']; % augmented array of derivatives of Tersoff potential in X-Y-Z directions

%**************************************************************************

function [fR,Dfrij] = fr(r,iMov,jBond)
% global x y z;
global A  B lamda1 lamda2 lamda3 alpha beta n c d h R D;
global DrijDxi DrijDyi DrijDzi

fR = A.*exp(-1.*lamda1.*r);
dfRdr = -1*lamda1*A*exp(-1*lamda1*r);

% DrijDxi= (x(iMov)-x(jBond))./r;
% DrijDyi= (y(iMov)-y(jBond))./r;
% DrijDzi= (z(iMov)-z(jBond))./r;

Dfrij(1,1)= dfRdr * DrijDxi;%force on atom i in X-direction
Dfrij(1,2)= dfRdr * DrijDyi;
Dfrij(1,3)= dfRdr * DrijDzi;
Dfrij(2,1)= -1*dfRdr * DrijDxi;%force on atom j in X-direction
Dfrij(2,2)= -1*dfRdr * DrijDyi;
Dfrij(2,3)= -1*dfRdr * DrijDzi;
%-------------------------------------------------------

function [fA,Dfaij] = fa(r,iMov,jBond)
global A  B lamda1 lamda2 lamda3 alpha beta n c d h R D;
global DrijDxi DrijDyi DrijDzi

fA = -1.*B.*exp(-1.*lamda2.*r);
dfAdr = -1.*lamda2.*-1.*B.*exp(-1.*lamda2.*r);

Dfaij(1,1)= dfAdr * DrijDxi;%force on atom i in X-direction
Dfaij(1,2)= dfAdr * DrijDyi;
Dfaij(1,3)= dfAdr * DrijDzi;
Dfaij(2,1)= -1*dfAdr * DrijDxi;%force on atom j in X-direction
Dfaij(2,2)= -1*dfAdr * DrijDyi;
Dfaij(2,3)= -1*dfAdr * DrijDzi;
%-------------------------------------------------------

function [fC,Dfcij] = fc(r)
global A  B lamda1 lamda2 lamda3 alpha beta n c d h R D;
global DrijDxi DrijDyi DrijDzi

if (r < (R-D))
	fC = 1;
	% end
elseif (r >= (R-D) & r <= (R+D))
	%                     fC(j,k,i)=0.5-0.5.*sin(pi./2.*(r(j,k,i)-R)./D);
	fC = 0.5+0.5.*cos(pi.*(r-(R-D))./(2.*D));
	% end
elseif (r > (R+D))
	fC=0;
end


if (r < (R-D))
	dfCdr = 0;
	% end
elseif (r >= (R-D) & r <= (R+D))
	dfCdr = -0.5.*(sin(((r-(R-D))./(2.*D)).*pi)).*(pi./(2.*D));
	%                     dfCdr(j,k,i)=-1/4.*cos(pi./2.*(R-r(j,k,i))./D).*pi./D;
	% end
elseif (r > (R+D))
	dfCdr=0;
end

Dfcij(1,1)= dfCdr * DrijDxi;
Dfcij(1,2)= dfCdr * DrijDyi;
Dfcij(1,3)= dfCdr * DrijDzi;
Dfcij(2,1)= -1*dfCdr * DrijDxi;
Dfcij(2,2)= -1*dfCdr * DrijDxi;
Dfcij(2,3)= -1*dfCdr * DrijDxi;
%--------------------------------------------------------------------------

function [g_ijk,Dg] = g(rij,rik,rjk)
global A  B lamda1 lamda2 lamda3 alpha beta n c d h R D;
global DrijDxi DrijDyi DrijDzi DrikDxi DrikDyi DrikDzi DrjkDxj DrjkDyj DrjkDzj


DrijDxj= -1*DrijDxi;
DrijDyj= -1*DrijDyi;
DrijDzj= -1*DrijDzi;

DrikDxk= -1*DrikDxi;
DrikDyk= -1*DrikDyi;
DrikDzk= -1*DrikDzi;

DrjkDxk= -1*DrjkDxj;
DrjkDyk= -1*DrjkDyj;
DrjkDzk= -1*DrjkDzj;


cosTh = (rij^2 + rik^2 - rjk^2)/(2 * rij * rik);

th=acos(cosTh);
sinTh = sin(th);

g_ijk = 1+ c^2/d^2 - c^2/(d^2+(h-cosTh)^2);

dgdth=2*c^2*((h-cosTh)*sinTh)/((d^2+(h-cosTh)^2)^2);

rij2=rij^2;
rik2=rik^2;
rjk2=rjk^2;

dthdrij = 1/2*(rij2-rik2+rjk2)/((rij2*rik)*(-1*sinTh));
dthdrik = 1/2*(rik2-rij2+rjk2)/((rik2*rij)*(-1*sinTh));
dthdrjk = -1*rjk/((rij*rik)* (-1*sinTh));


dthdxi = dthdrij*DrijDxi + dthdrik*DrikDxi;% calculate components w.r.t. X-Y-Z coordinates of atoms i,j,k 
dthdyi = dthdrij*DrijDyi + dthdrik*DrikDyi;
dthdzi = dthdrij*DrijDzi + dthdrik*DrikDzi;

dthdxj = dthdrij*DrijDxj + dthdrjk*DrjkDxj;
dthdyj = dthdrij*DrijDyj + dthdrjk*DrjkDyj;
dthdzj = dthdrij*DrijDzj + dthdrjk*DrjkDzj;

dthdxk = dthdrik*DrikDxk + dthdrjk*DrjkDxk;
dthdyk = dthdrik*DrikDyk + dthdrjk*DrjkDyk;
dthdzk = dthdrik*DrikDzk + dthdrjk*DrjkDzk;

Dg(1,1) = dgdth*dthdxi;
Dg(1,2) = dgdth*dthdyi;
Dg(1,3) = dgdth*dthdzi;

Dg(2,1) = dgdth*dthdxj;
Dg(2,2) = dgdth*dthdyj;
Dg(2,3) = dgdth*dthdzj;

Dg(3,1) = dgdth*dthdxk;
Dg(3,2) = dgdth*dthdyk;
Dg(3,3) = dgdth*dthdzk;
%-------------------------------------------------------------------

function [zetaij,Dzetaij,Dzetak] = zeta(k3Body,rij,rik,rjk,zetaij,Dzetaij,Dzetak)

global A  B lamda1 lamda2 lamda3 alpha beta n c d h R D;
global DrijDxi DrijDyi DrijDzi DrikDxi DrikDyi DrikDzi DrjkDxj DrjkDyj DrjkDzj

% zetaij=0;
% Dzetaij=zeros(2,3);
% 		
% for k=2:list(iMov,1)+1
% 	k3Body=list(iMov,k);
% 	if((k3Body == iMov) | (k3Body == jBond))
% 		continue
% 	end
% 	rik= sqrt( (x(iMov)-x(k3Body)).^2 + (y(iMov)-y(k3Body)).^2 + (z(iMov)-z(k3Body)).^2 );
% 	rjk= sqrt( (x(jBond)-x(k3Body)).^2 + (y(jBond)-y(k3Body)).^2 + (z(jBond)-z(k3Body)).^2 );
	
	
% 	DrikDxi= (x(iMov)-x(k3Body))./rik;
% 	DrikDyi= (y(iMov)-y(k3Body))./rik;
% 	DrikDzi= (z(iMov)-z(k3Body))./rik;
% 	
% 	DrjkDxj= (x(jBond)-x(k3Body))./rjk;%????????????????????????????
% 	DrjkDyj= (y(jBond)-y(k3Body))./rjk;
% 	DrjkDzj= (z(jBond)-z(k3Body))./rjk;
	
	
	[fC,Dfcik] = fc(rik);
	
	[g_ijk,Dg] = g(rij,rik,rjk);
	
	zetaij = zetaij + fC*g_ijk;
	
	Dzetaij(1,1) = Dzetaij(1,1) + (fC*Dg(1,1) + Dfcik(1,1)*g_ijk);%derivative w.r.r. xi
	Dzetaij(1,2) = Dzetaij(1,2) + (fC*Dg(1,2) + Dfcik(1,2)*g_ijk);
	Dzetaij(1,3) = Dzetaij(1,3) + (fC*Dg(1,3) + Dfcik(1,3)*g_ijk);
	
	Dzetaij(2,1) = Dzetaij(2,1) + fC*Dg(2,1);%derivative w.r.r. xj
	Dzetaij(2,2) = Dzetaij(2,2) + fC*Dg(2,2);
	Dzetaij(2,3) = Dzetaij(2,3) + fC*Dg(2,3);
	
	Dzetak(k3Body,1)= fC*Dg(3,1) + g_ijk*Dfcik(2,1);%derivative w.r.r. xk
	Dzetak(k3Body,2)= fC*Dg(3,2) + g_ijk*Dfcik(2,2);
	Dzetak(k3Body,3)= fC*Dg(3,3) + g_ijk*Dfcik(2,3);
	
	% 				end
% end
%--------------------------------------------------------------------------

function [bij,Dbij,Dbk] = b(iMov,jBond,zetaij,Dzetaij,Dzetak,list,Dbij,Dbk)

global A  B lamda1 lamda2 lamda3 alpha beta eta c d h R D;

bij=(1+ (beta^eta * zetaij^eta))^(-1/(2*eta));

if(zetaij == 0)
	zetaij = 1.0e-10; % to avoid infinity if zeta is 0 then in the calculation of derivative of the bij index of beta becomes -ve, therefore it becomes infinity
end
		
% 			if(zeta ~= 0)
Dbij(1,1) = -1/2*((1+(beta^eta)*(zetaij^eta))^(-1/2*(1+2*eta)/eta))*(beta^eta)*zetaij^(eta-1) .* Dzetaij(1,1);
Dbij(1,2) = -1/2*((1+(beta^eta)*(zetaij^eta))^(-1/2*(1+2*eta)/eta))*(beta^eta)*zetaij^(eta-1) .* Dzetaij(1,2);
Dbij(1,3) = -1/2*((1+(beta^eta)*(zetaij^eta))^(-1/2*(1+2*eta)/eta))*(beta^eta)*zetaij^(eta-1) .* Dzetaij(1,3);

Dbij(2,1) = -1/2*((1+(beta^eta)*(zetaij^eta))^(-1/2*(1+2*eta)/eta))*(beta^eta)*zetaij^(eta-1) .* Dzetaij(2,1);
Dbij(2,2) = -1/2*((1+(beta^eta)*(zetaij^eta))^(-1/2*(1+2*eta)/eta))*(beta^eta)*zetaij^(eta-1) .* Dzetaij(2,2);
Dbij(2,3) = -1/2*((1+(beta^eta)*(zetaij^eta))^(-1/2*(1+2*eta)/eta))*(beta^eta)*zetaij^(eta-1) .* Dzetaij(2,3);

for k=2:list(iMov,1)+1
	k3Body=list(iMov,k);
	if((k3Body == iMov) | (k3Body == jBond))
		continue;
	end
	Dbk(k3Body,1) = -1/2*((1+(beta^eta)*(zetaij^eta))^(-1/2*(1+2*eta)/eta))*(beta^eta)*zetaij^(eta-1) .* Dzetak(k3Body,1);
	Dbk(k3Body,2) = -1/2*((1+(beta^eta)*(zetaij^eta))^(-1/2*(1+2*eta)/eta))*(beta^eta)*zetaij^(eta-1) .* Dzetak(k3Body,2);
	Dbk(k3Body,3) = -1/2*((1+(beta^eta)*(zetaij^eta))^(-1/2*(1+2*eta)/eta))*(beta^eta)*zetaij^(eta-1) .* Dzetak(k3Body,3);
	% 				end
end
