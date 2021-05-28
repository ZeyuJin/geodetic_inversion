function [ux,uy,uz]=calc_okada(HF,U,x,y,nu,delta,d,len,W,fault_type,strike,tp);
%[ux,uy,uz]=calc_okada(HF,U,x,y,nu,delta,d,len,W,fault_type,strike,tp);
% HF is a scaling factor for horizontal displacements (e.g., anisotropy)
% U is slip
% x,y are the observation points
% d is depth (positive down)
% nu is Poisson ratio
% delta is dip angle
% len,W are the fault length and width, resp.
% fault_type is 1 2 3 for strike, dip, and opening
% tp is the optional topo vector (size(x))
%

cosd  = cos(delta);
sind  = sin(delta);

%% define parameters with respect to the BOTTOM of the fault (assume
%% that the input ones correspond to the TOP)

d=d+W*sind;
x=x-W*cosd*cos(strike); 
y=y+W*cosd*sin(strike);

d=d+tp;
strike = -strike+pi/2;
coss  = cos(strike);
sins  = sin(strike);
rot = [coss -sins ; sins coss];
rotx =  x*coss+y*sins;
roty = -x*sins+y*coss;

%%%%% Okada fault model for dislocation in an elastic half-space.
%%%%% based on BSSA Vol. 95 p.1135-45, 1985

L     = len/2;
Const = -U/(2*pi);

p = roty*cosd + d*sind;	%a matrix eqn. (30)
q = roty*sind - d*cosd;	%a matrix eqn. (30)
a = 1-2*nu;		% mu/(lambda+mu) = 1-2*poisson's ratio

parvec = [a, delta, fault_type];

[f1a,f2a,f3a] = fBi(rotx+L, p  , parvec, p, q);
[f1b,f2b,f3b] = fBi(rotx+L, p-W, parvec, p, q);
[f1c,f2c,f3c] = fBi(rotx-L, p  , parvec, p, q);
[f1d,f2d,f3d] = fBi(rotx-L, p-W, parvec, p, q);

%%%%% Displacement eqns. (25-27)

uxj = Const * (f1a - f1b - f1c + f1d);
uyj = Const * (f2a - f2b - f2c + f2d);
uz = Const * (f3a - f3b - f3c + f3d);

% rotate horizontals back to the orig. coordinate system
ux=HF*(-uyj*sins+uxj*coss);  
uy=HF*(uxj*sins+uyj*coss);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%END OF THE calc_okada.m%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f1,f2,f3] = fBi(sig, eta, parvec, p, q);

%%%%% Pull out parameters form parvec    %%%%%
%%%%% parvec = [a, delta, fault_type] %%%%%

a          = parvec(1);
delta      = parvec(2);
fault_type = parvec(3);

%%%%% Some simplifications %%%%%

epsn  = 1.0e-15;
cosd  = cos(delta);
sind  = sin(delta);
tand  = tan(delta);
cosd2 = cos(delta)^2;
sind2 = sin(delta)^2;
cssnd = cos(delta)*sin(delta);

R     = sqrt(sig.^2 + eta.^2 + q.^2);
X     = sqrt(sig.^2 + q.^2);		
ytil  = eta*cosd + q*sind;		
dtil  = eta*sind - q*cosd;		

Rdtil = R + dtil;			
Rsig  = R + sig;			
Reta  = R + eta;			
RX    = R + X;
			
lnRdtil = log(Rdtil);		
lnReta  = log(Reta);	
lnReta0 = -log(R-eta);	
ORRsig  = 1 ./ (R .* Rsig);			
OReta   = 1 ./ Reta;			
ORReta  = 1 ./ (R .* Reta);

indfix  = find(abs(Reta) < epsn);		%check for bad values
if (~isempty(indfix))
  lnReta(indfix) = lnReta0(indfix);		
  OReta(indfix)  = 0 * indfix;
  ORReta(indfix) = 0 * indfix;
end

indfix  = find(abs(Rsig) < epsn);		
if (~isempty(indfix))
  ORsig(indfix)  = 0 * indfix;
  ORRsig(indfix) = 0 * indfix;
end

%%%%% theta term with q = 0 fix

%betatop=sig.*eta;
%betabottom=q.*R;
%nz=find(betabottom~=0);
%theta=pi/2*sign(betatop);
%theta(nz)=atan(betatop(nz)./betabottom(nz));
%theta(find(abs(q) > epsn))  = atan((sig.*eta)./(q.*R));		
indfix = find(abs(q) <= epsn);
if (~isempty(indfix))
  theta(indfix) = 0 * indfix;
  indfix = find(abs(q) > epsn);
  theta(indfix) = atan((sig(indfix).*eta(indfix))./(q(indfix).*R(indfix)));
  theta=reshape(theta,size(sig));

 else
  theta = atan((sig.*eta)./(q.*R));
end

%%%%% The I_12345 factors %%%%%

if abs(cosd) < epsn
%%%%% cosd = 0 fix [eqn. (29)]
  I5 = -a   .* sig .* sind ./ Rdtil;
  I4 = -a   .* q ./ Rdtil;
  I3 =  a/2 .* (eta ./ Rdtil + (ytil .* q) ./ (Rdtil.^2) - lnReta );
  I2 = -a   .* lnReta - I3;
  I1 = -a/2 .* (sig .* q) ./ (Rdtil.^2);
else
%%%%% default [eqn. (28)]
%   sigtemp(indfix) = epsn;
%  I5 = a * 2 ./ cosd .* atan2( (eta.*(X+q.*cosd) + X.*RX.*sind),...
%                             (sig.*RX.*cosd) );
%%  indfix = find(abs(sig)<epsn);
   I5 = a * 2 ./ cosd .* ...
     atan( (eta.*(X+q.*cosd) + X.*RX.*sind)./(sig.*RX.*cosd) ); 
  indfix = find(abs(sig)<epsn);
  if (~isempty(indfix))
   I5(indfix) = 0 * indfix;
  end
  I4 = a ./ cosd .* (lnRdtil - sind .* lnReta);
  I3 = a * (1 ./ cosd .* ytil ./ Rdtil - lnReta) + tand .* I4;
  I2 = -a .* lnReta - I3;
  I1 = -a ./ cosd .* sig ./ Rdtil - tand .* I5;
end

%%%%% The fault specific parameters %%%%%

if fault_type == 1;  		%%%%% Strike Slip [eqn. (25)]
%2+2
%size(theta)

  f1 = (sig .* q)  .* ORReta + theta + I1 .* sind;
  f2 = (ytil .* q) .* ORReta + (q .* cosd) .* OReta + I2 .* sind;
  f3 = (dtil .* q) .* ORReta + (q .* sind) .* OReta + I4 .* sind;

elseif fault_type == 2;  	%%%%% Dip Slip  [eqn. (26)]

  f1 = q./R - I3 .* cssnd;
  f2 = (ytil .* q) .* ORRsig + cosd .* theta - I1 .* cssnd;
  f3 = (dtil .* q) .* ORRsig + sind .* theta - I5 .* cssnd;

else fault_type == 3;		%%%%% Tensile [eqn. (27)]

  f1 = q.^2 .* ORReta - I3 .* sind2;
  f2 = (-dtil .* q) .* ORRsig ...
     - sind .* ((sig .* q) .* ORReta - theta) - I1 .* sind2;
  f3 = (ytil .* q) .* ORRsig ...
     + cosd .* ((sig .* q) .* ORReta - theta) - I5 .* sind2;

end

