%% Set initial condition, target, blocks, parameters
clear; close all; clc
load batteryParams.mat
load info0.mat
params = struct('m', 1348 + 4*80 + 4*250*0.895,... % mass of eVTOL 2563kg 
                'rou', 1.225,... % air density kg/m^3
                'Cd', 1,...
                'Fx', 1.5,...   %2.11 for ehang184? Paper DE2 f is CD*Fx
                'Fh', 4,...
                'g', 9.8,...
                'A', 24.6,...
                'vtip', 120,...
                'eWidth', 8,...
                'eHeight', 3,...
                'mBat',0.895,...
                'cBat',1380,...
                'nBat',4*250,...
                'ParaOCV',ParaOCV,... % [6,1]
                'ParaRo',ParaRo,... % [14,1]
                'ParaRp',ParaRp,... % [14,1]
                'ParaCp',ParaCp,... % [14,1]
                'ParadUdT',ParadUdT,... % [7,1]
                'ObsX',10000,...
                'ObsH',0,...
                'ObsWidth',20000-2000*2,...
                'ObsHeight',600-40,...
                'sigma',0.1,...
                'Cd0',0.01);

% states of the State Function:
% [vx; vh; x; h; SOC; Up; Tbat; I]
p1 = 25; p2 = 76; p3 = 20;

targetPose1 = [50;0;inf;310;inf;inf;inf;inf];
p = zeros(7,1);
Ts = zeros(7,1);
rectangle('Position',[2000,0,20000-4000,280]); hold on; axis equal
%% Stage 1 Take off
initialPose = [0;0;0;0;1;0;25+273.15;0];
targetPose = [0;0;0;15;inf;inf;inf;inf];
u0 = [params.m*params.g;0];
p(1) = 7; Ts(1) = 1;
[msobj,simdata] = eVTOLmsobj(p(1),Ts(1),initialPose,targetPose,u0,params,'eVTOLcostTakeoff');
tic;[~,~,info{1}] = nlmpcmove(msobj,initialPose,u0,simdata);t=toc;
fprintf('Stage1>>> Calculation Time = %s; Objective cost = %s; ExitFlag = %s; Iterations = %s\n',...
    num2str(t),num2str(info{1}.Cost),num2str(info{1}.ExitFlag),num2str(info{1}.Iterations));
plot(info{1}.Xopt(1:end-2,3),info{1}.Xopt(1:end-2,4),'-o'); hold on; axis equal

%% Stage 2 Hover
initialPose = info{1}.Xopt(p(1)-1,:);
targetPose = [0;0;0;15;inf;inf;inf;inf];
u0 = info{1}.MVopt(p(1)-1,:);
p(2) = 30; Ts(2) = 2;
[msobj,simdata] = eVTOLmsobj(p(2),Ts(2),initialPose,targetPose,u0,params,'eVTOLcostHover1');
tic;[~,~,info{2}] = nlmpcmove(msobj,initialPose,u0,simdata);t=toc;
fprintf('Stage2>>> Calculation Time = %s; Objective cost = %s; ExitFlag = %s; Iterations = %s\n',...
    num2str(t),num2str(info{2}.Cost),num2str(info{2}.ExitFlag),num2str(info{2}.Iterations));
plot(info{2}.Xopt(1:end-2,3),info{2}.Xopt(1:end-2,4),'-o'); hold on; axis equal

%% Stage 3 Climb
initialPose = info{2}.Xopt(p(2)-1,:);
targetPose = [44.2;0;2512;300;inf;inf;inf;inf];
u0 = info{2}.MVopt(p(2)-1,:);
p(3) = 22; Ts(3) = 3;
% [msobj,simdata] = eVTOLmsobj(p(3),Ts(3),initialPose,targetPose,u0,params,'eVTOLcostClimb','eVTOLcostClimbLast');
[msobj,simdata] = eVTOLmsobjClimb(p(3),Ts(3),initialPose,targetPose,u0,params,'eVTOLcostClimb');
tic;[~,~,info{3}] = nlmpcmove(msobj,initialPose,u0,simdata);t=toc;
fprintf('Stage3>>> Calculation Time = %s; Objective cost = %s; ExitFlag = %s; Iterations = %s\n',...
    num2str(t),num2str(info{3}.Cost),num2str(info{3}.ExitFlag),num2str(info{3}.Iterations));
plot(info{3}.Xopt(1:end,3),info{3}.Xopt(1:end,4),'-o'); hold on; axis equal

%% Stage 4 Cruise
initialPose = info{3}.Xopt(p(3)+1,:);
targetPose = [44.2;0;17098;300;inf;inf;inf;inf];
u0 = info{3}.MVopt(p(3)-1,:);
p(4) = 33; Ts(4) = 10;
[msobj,simdata] = eVTOLmsobj(p(4),Ts(4),initialPose,targetPose,u0,params,'eVTOLcostCruise');
tic;[~,~,info{4}] = nlmpcmove(msobj,initialPose,u0,simdata);t=toc;
fprintf('Stage4>>> Calculation Time = %s; Objective cost = %s; ExitFlag = %s; Iterations = %s\n',...
    num2str(t),num2str(info{4}.Cost),num2str(info{4}.ExitFlag),num2str(info{4}.Iterations));
plot(info{4}.Xopt(1:end-2,3),info{4}.Xopt(1:end-2,4),'-o'); hold on; axis equal
% Cal_Energy(info{4},Ts(4),params)
% Cal_Energy(info0{4},1,params)

%% Stage 5 Descent
% initialPose = info{4}.Xopt(p(4)-1,:);
initialPose = info{4}.Xopt(p(4)+1,:);
targetPose = [0;0;20000;15;inf;inf;inf;inf];
u0 = info{4}.MVopt(p(4),:);
p(5) = 26; Ts(5) = 3;
[msobj,simdata] = eVTOLmsobjDescent(p(5),Ts(5),initialPose,targetPose,u0,params,'eVTOLcostDescent');
tic;[~,~,info{5}] = nlmpcmove(msobj,initialPose,u0,simdata);t=toc;
fprintf('Stage5>>> Calculation Time = %s; Objective cost = %s; ExitFlag = %s; Iterations = %s\n',...
    num2str(t),num2str(info{5}.Cost),num2str(info{5}.ExitFlag),num2str(info{5}.Iterations));
plot(info{5}.Xopt(1:end,3),info{5}.Xopt(1:end,4),'-o'); hold on; axis equal
% Cal_Energy(info{5},Ts(5),params)
% Cal_Energy(info0{5},1,params)
% Cal_aging(info{5},Ts(5))
% Cal_aging(info0{5},1)
% save('Descent.mat')
%% Stage 6 Hover
initialPose = info{5}.Xopt(p(5)+1,:);
targetPose = [0;0;20000;15;inf;inf;inf;inf];
u0 = info{5}.MVopt(p(5)-1,:);
p(6) = 30; Ts(6) = 2;
[msobj,simdata] = eVTOLmsobj(p(6),Ts(6),initialPose,targetPose,u0,params,'eVTOLcostHover2');
tic;[~,~,info{6}] = nlmpcmove(msobj,initialPose,u0,simdata);t=toc;
fprintf('Stage6>>> Calculation Time = %s; Objective cost = %s; ExitFlag = %s; Iterations = %s\n',...
    num2str(t),num2str(info{6}.Cost),num2str(info{6}.ExitFlag),num2str(info{6}.Iterations));
plot(info{6}.Xopt(1:end-2,3),info{6}.Xopt(1:end-2,4),'-o'); hold on; axis equal
% Cal_Energy(info{6},Ts(6),params)
% Cal_Energy(info0{6},1,params)
% Cal_aging(info{6},Ts(6))
% Cal_aging(info0{6},1)
% save('Hover2.mat')
%% Stage 7 Landing
initialPose = info{6}.Xopt(p(6)-1,:);
targetPose = [0;0;20000;0;inf;inf;inf;inf];
u0 = info{6}.MVopt(p(6)-1,:);
p(7) = 8; Ts(7) = 1;
[msobj,simdata] = eVTOLmsobj(p(7),Ts(7),initialPose,targetPose,u0,params,'eVTOLcostLanding');
tic;[~,~,info{7}] = nlmpcmove(msobj,initialPose,u0,simdata);t=toc;
fprintf('Stage1>>> Calculation Time = %s; Objective cost = %s; ExitFlag = %s; Iterations = %s\n',...
    num2str(t),num2str(info{7}.Cost),num2str(info{7}.ExitFlag),num2str(info{7}.Iterations));
plot(info{7}.Xopt(1:end-2,3),info{7}.Xopt(1:end-2,4),'-o'); hold on; axis equal
% Cal_Energy(info{7},Ts(7),params)
% Cal_Energy(info0{7},1,params)
% Cal_aging(info{7},Ts(7))
% Cal_aging(info0{7},1)
% save('Finish.mat')
%%
[Xopt,MVopt,Topt,Popt] = eVTOLplot(initialPose, targetPose1, params, info);
% save('info.mat',"info") 
% save('Xopt.mat',"Xopt") 
% save('MVopt.mat',"MVopt") 
% save('Topt.mat',"Topt") 
% save('Popt.mat',"Popt") 
% save('NP.mat');