function cost =  eVTOLcostClimb(stage,x,u,dmv,p)

% Copyright 2020-2023 The MathWorks, Inc.

%#codegen
% Wmv = eye(2);
Wdmv = diag([1,10^7]);
% CostMv = u'*Wmv*u; 
CostDmv = dmv'*Wdmv*dmv;
% cost = CostDmv+CostMv;


rou = p(2);
Area = p(7);
nBat = p(65);
vtip = p(66);
sigma = p(67);
Cd0 = p(68);
T = u(1);
theta = u(2);

v = sqrt(x(1).^2 + x(2).^2);
gamma = atan2(x(2),x(1));
alpha = gamma + theta;
vi0 = sqrt(T/2/rou/Area);
funVi = @(vi) vi-vi0^2/sqrt(v^2+2*v*sin(alpha)*vi+vi^2);
vi = max(fsolve(funVi, 100, optimset('Display', 'off')));
Pinduced = T*(vi*1.15+v*sin(alpha));
Pprofile = (1+4.7*(x(1)/vtip).^2)*rou*Area*vtip^3*sigma/8*Cd0;
P = (Pinduced+Pprofile)/nBat; 
% cost = 0.01*P;

cost = P+0.5*CostDmv;