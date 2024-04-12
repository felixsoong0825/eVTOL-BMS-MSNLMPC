function dxdt = eVTOLstateFcn(x,u,p)
vx = x(1);  
vh = x(2);
SOC = x(5);
Up = x(6);
TempK = x(7);
TempC = x(7) - 273.15;
I_pre = x(8);

T = u(1);
theta = u(2);

m = p(1);
rou = p(2);
Cd = p(3);
Fx = p(4);
Fh = p(5);
g = p(6);
Area = p(7);
ParaOCV = p(8:8+5);
ParaRo = p(14:14+13);
ParaRp = p(28:28+13);
ParaCp = p(42:42+13);
ParadUdT = p(56:56+6);
mBat = p(63);
cBat = p(64);
nBat = p(65);
vtip = p(66);
sigma = p(67);
Cd0 = p(68);
Ts = p(69);

v = sqrt(x(1).^2 + x(2).^2);
gamma = atan2(x(2),x(1));
alpha = gamma + theta;
vi0 = sqrt(T/2/rou/Area);
funVi = @(vi) vi-vi0^2/sqrt(v^2+2*v*sin(alpha)*vi+vi^2);
vi = max(fsolve(funVi, 100, optimset('Display', 'off')));
Pinduced = T*(vi*1.15+v*sin(alpha));
Pprofile = (1+4.7*(x(1)/vtip).^2)*rou*Area*vtip^3*sigma/8*Cd0;
P = (Pinduced+Pprofile)/0.9/nBat; 

Uoc = [SOC^5 SOC^4 SOC^3 SOC^2 SOC 1]*ParaOCV;
% funRo = @(Crate) 0.001*[1 TempC TempC.^2 TempC.*Crate SOC SOC.^2 SOC.^3 SOC.^4 SOC.*TempC Crate Crate.^2 Crate.^3 Crate.^4 Crate.*SOC]*ParaRo;
% funCrate = @(Crate) Crate*50-(Uoc-Up-((Up-Uoc)^2-4*funRo(Crate)*P)^0.5)/2/funRo(Crate);
% Crate = fsolve(funCrate, 0, optimset('Display', 'off'));
% Crate = I_pre/50;
Crate = 3.75;
Ro = 0.001*[1 TempC TempC.^2 TempC.*Crate SOC SOC.^2 SOC.^3 SOC.^4 SOC.*TempC Crate Crate.^2 Crate.^3 Crate.^4 Crate.*SOC]*ParaRo;
Rp = 0.001*[1 TempC TempC.^2 TempC.*Crate SOC SOC.^2 SOC.^3 SOC.^4 SOC.*TempC Crate Crate.^2 Crate.^3 Crate.^4 Crate.*SOC]*ParaRp;
Cp =  [1 TempC TempC.^2 TempC.*Crate SOC SOC.^2 SOC.^3 SOC.^4 SOC.*TempC Crate Crate.^2 Crate.^3 Crate.^4 Crate.*SOC]*ParaCp;
dUdT = 0.001 * [1 SOC SOC^2 SOC^3 SOC^4 SOC^5 SOC^6]*ParadUdT;
% I = min(roots([Ro,Up-Uoc,P])); % I^2*Ro + I*(Up-Uoc) + P = 0
I = (Uoc-Up-((Up-Uoc)^2-4*Ro*P)^0.5)/2/Ro;
% I = Crate*50;

dxdt = zeros(5,1);
dxdt(1) = -0.5/m*rou*Cd*Fx*vx*abs(vx)+1/m*T*sin(theta);
dxdt(2) = -0.5/m*rou*Cd*Fh*vh*abs(vh)+1/m*T*cos(theta)-g;
% dxdt(1) = -0.5/m*rou*Cd*Fx*vx^2+1/m*T*sin(theta);
% dxdt(2) = -0.5/m*rou*Cd*Fh*vh^2+1/m*T*cos(theta)-g;
dxdt(3) = vx;
dxdt(4) = vh;
dxdt(5) = -I/50/3600; % I 'discharge' = '+'
dxdt(6) = I/Cp-Up/Rp/Cp; 
dxdt(7) = (I^2*(Ro+Rp)+I*TempK*dUdT-(1.478*(TempC-25)))/mBat/cBat;
dxdt(8) = (I-I_pre)/Ts;