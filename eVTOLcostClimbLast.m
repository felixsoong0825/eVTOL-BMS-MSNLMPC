function cost =  eVTOLcostClimbLast(stage,x,u,dmv,p)
% Cost function of the path planner of a truck and trailer system

% Copyright 2020-2023 The MathWorks, Inc.

%#codegen
% Wmv = eye(2);
Wdmv = eye(2);
% Wx = diag([0,1,-1,1]);
% CostMv = u'*Wmv*u; 
CostDmv = dmv'*Wdmv*dmv;
% CostX = (x(1:4)-[0;0;0;310])'*Wx*(x(1:4)-[0;0;0;310]);
cost = CostDmv;

