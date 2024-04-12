function [msobj,simdata] =  eVTOLmsobjClimb(p,Ts,initialPose,targetPose,u0,params,costfun,costfunLast)
msobj = nlmpcMultistage(p,8,2);
msobj.Ts = Ts;
msobj.UseMVRate = true;

% msobj.MV(1).Min = 0; 
% msobj.MV(1).Max = 30000;     
msobj.MV(1).RateMin = -msobj.Ts*2000; msobj.MV(1).RateMax = msobj.Ts*2000;
msobj.MV(2).Min = -pi/7; msobj.MV(2).Max =  pi/7; 
msobj.MV(2).RateMin = -5*msobj.Ts*pi/180; msobj.MV(2).RateMax = 5*msobj.Ts*pi/180;

msobj.States(1).Min = 0; msobj.States(1).Max =  44.2;
msobj.States(2).Min = 0; msobj.States(2).Max =  15;
% msobj.States(3).Min = 0; msobj.States(3).Max = 30000;
% msobj.States(4).Min = 0; msobj.States(4).Max = 400;
% msobj.States(5).Min = 0; msobj.States(5).Max = 1;
% msobj.States(6).Min = 0; msobj.States(6).Max = 4;
msobj.States(7).Min = 0; msobj.States(7).Max = 323.15;
msobj.States(8).Min = 0; msobj.States(8).Max = 250;

msobj.Model.TerminalState = targetPose;
msobj.Model.ParameterLength = 68+1;
msobj.Model.StateFcn = 'eVTOLstateFcn';

for ct=1:p
    msobj.Stages(ct).CostFcn = costfun;
    msobj.Stages(ct).ParameterLength = 68+8;
    msobj.Stages(ct).IneqConFcn = "eVTOLineqConFcn";
end

if nargin > 7
    for ct=p
        msobj.Stages(ct).CostFcn = costfunLast;
    end
end

simdata = getSimulationData(msobj,'TerminalState');
simdata.StateFcnParameter = [params.m;params.rou;params.Cd;params.Fx;params.Fh;params.g;params.A;params.ParaOCV;params.ParaRo;params.ParaRp;params.ParaCp;params.ParadUdT;params.mBat;params.cBat;params.nBat;params.vtip;params.sigma;params.Cd0;Ts];
simdata.StageParameter = repmat([params.m;params.rou;params.Cd;params.Fx;params.Fh;params.g;params.A;params.ParaOCV;params.ParaRo;params.ParaRp;params.ParaCp;params.ParadUdT;params.mBat;params.cBat;params.nBat;params.vtip;params.sigma;params.Cd0;params.eWidth;params.eHeight;params.A;params.rou;params.ObsX;params.ObsH;params.ObsWidth;params.ObsHeight],p,1);
simdata.TerminalState = targetPose;

validateFcns(msobj,[2;3;0.5;0.4;1;0;300;50],[0;0.2],simdata);
fprintf('Automated Flying Planner is running...\n');

cineq = eVTOLineqConFcn(1,initialPose,u0,0,...
    [params.m;params.rou;params.Cd;params.Fx;params.Fh;params.g;params.A;params.ParaOCV;params.ParaRo;params.ParaRp;params.ParaCp;params.ParadUdT;params.mBat;params.cBat;params.nBat;params.vtip;params.sigma;params.Cd0;params.eWidth;params.eHeight;params.A;params.rou;params.ObsX;params.ObsH;params.ObsWidth;params.ObsHeight]);
if any(cineq>0)
    fprintf('Initial pose is not valid.\n');
    return
end