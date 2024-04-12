function [Xopt,MVopt,Topt,Popt] = eVTOLplot(initialPose, targetPose, params, info)
%EVTOLPLOT 此处显示有关此函数的摘要
%   此处显示详细说明

Xopt = [info{1}.Xopt(1:end-3,:);info{2}.Xopt(1:end-3,:);info{3}.Xopt(1:end-1,:);info{4}.Xopt(1:end-1,:);info{5}.Xopt(1:end-1,:);info{6}.Xopt(1:end-3,:);info{7}.Xopt];
MVopt = [info{1}.MVopt(1:end-3,:);info{2}.MVopt(1:end-3,:);info{3}.MVopt(1:end-1,:);info{4}.MVopt(1:end-1,:);info{5}.MVopt(1:end-1,:);info{6}.MVopt(1:end-3,:);info{7}.MVopt];

for i = 1:6
    topt{i} = info{i}.Topt(1:end-2);
end
topt{3} = info{3}.Topt;
topt{4} = info{4}.Topt;
topt{5} = info{5}.Topt;
topt{7} = info{7}.Topt;

Topt = topt{1};
for i = 2:7
   Topt = [Topt(1:end-1);Topt(end)+topt{i}];
end

ct = length(Xopt);
%% set
close all
f = figure('NumberTitle','off');
f.Name = 'EVTOL Trajectory with Initial and Target Positions';
xlabel('x'); ylabel('y'); hold on;
xlim([17500 21000])
axis equal
%%
ObsBox = collisionBox(params.ObsWidth,params.ObsHeight,0);
ObsBox.Pose(1:2,4) = [params.ObsX params.ObsH];
[~,Obsshow] = show(ObsBox);
Obsshow.FaceColor = 'y';
%%
% plot(XY0(:,1),XY0(:,2), '-o', MarkerSize=5)
plot(Xopt(:,3),Xopt(:,4), '-o', MarkerSize=5)
%%
Area = params.A;
rou = params.rou;
nBat = params.nBat;
vtip = params.vtip;
sigma = params.sigma;
Cd0 = params.Cd0;
ParaRo = params.ParaRo;
ParaRp = params.ParaRp;
ParaCp = params.ParaCp;
ParaOCV = params.ParaOCV;
vx = Xopt(:,1); vh = Xopt(:,2);
x = Xopt(:,3); h = Xopt(:,4);
SOC = Xopt(:,5); Up = Xopt(:,6); TempK = Xopt(:,7); TempC = TempK-273.15;
T = MVopt(:,1); theta = MVopt(:,2);
v = sqrt(vx.^2 + vh.^2);
gamma = atan2(vh,vx);
alpha = gamma+theta;
vi0 = sqrt(T/2/rou/Area);
vi = zeros(ct,1);
for i=1:ct
    fun = @(vi) vi-vi0(i)^2/sqrt(v(i)^2+2*v(i)*sin(alpha(i))*vi+vi^2);
    % fun = @(vi) vi^4+2*v(i)*sin(alpha(i))*vi^3+v(i)^2*vi^2-vi0(i)^4;
    vi(i) = max(fsolve(fun, 10, optimset('Display', 'off')));
end
Pinduced = 1.15*T.*vi./nBat;
Pclimb = T.*v.*sin(alpha)./nBat;
Pprofile = 2*ones(ct,1)*(rou*Area*vtip^3*sigma/8*Cd0)./nBat;
Power = Pinduced+Pprofile+Pclimb; 
Uoc = [SOC.^5 SOC.^4 SOC.^3 SOC.^2 SOC ones(ct,1)]*ParaOCV;
Crate = zeros(ct,1);
for i=1:ct    
    funRo = @(Crate) 0.001*[1 TempC(i) TempC(i).^2 TempC(i).*Crate SOC(i) SOC(i).^2 SOC(i).^3 SOC(i).^4 SOC(i).*TempC(i) Crate Crate.^2 Crate.^3 Crate.^4 Crate.*SOC(i)]*ParaRo;
    funCrate = @(Crate) Crate*50-(Uoc(i)-Up(i)-((Up(i)-Uoc(i))^2-4*funRo(Crate)*Power(i))^0.5)/2/funRo(Crate);
    Crate(i) = fsolve(funCrate, 1, optimset('Display', 'off'));
end
% Crate = 3.76*ones(ct,1);
Ro = 0.001*[ones(ct,1) TempC TempC.^2 TempC.*Crate SOC SOC.^2 SOC.^3 SOC.^4 SOC.*TempC Crate Crate.^2 Crate.^3 Crate.^4 Crate.*SOC]*ParaRo;
Rp = 0.001*[ones(ct,1) TempC TempC.^2 TempC.*Crate SOC SOC.^2 SOC.^3 SOC.^4 SOC.*TempC Crate Crate.^2 Crate.^3 Crate.^4 Crate.*SOC]*ParaRp;
Cp = [ones(ct,1) TempC TempC.^2 TempC.*Crate SOC SOC.^2 SOC.^3 SOC.^4 SOC.*TempC Crate Crate.^2 Crate.^3 Crate.^4 Crate.*SOC]*ParaCp;

% I = (Uoc-Up-((Up-Uoc).^2-4.*Ro.*power).^0.5)/2./Ro;
I = 50*Crate;
U = Uoc - I.*Ro - Up;
%%
for i = 1:ct
    plot(initialPose(3), initialPose(4),'rx'); 
    plot(targetPose(3), targetPose(4),'gx'); 
    eVTOLbox = collisionBox(params.eWidth, params.eHeight,0);
    eVTOLbox.Pose(1:3,1:3) = eul2rotm([-(theta(i)) 0 0]);
    eVTOLbox.Pose(1:2,4) = [x(i);h(i)];
    show(eVTOLbox)
    % pause(0.001)
end
%%
f = figure('NumberTitle','off');
f.Name = 'EVTOL State';
subplot(3,4,1); hold on; box on
plot(Topt,x, '-o', MarkerSize=3)
plot(Topt,h, '-o', MarkerSize=3)
xlabel('t'); ylabel('d');
legend('x','h')

subplot(3,4,2); hold on; box on
plot(Topt,vx)
plot(Topt,vh)
plot(Topt,vi)
plot(Topt,v)
xlabel('t'); ylabel('v');
legend('vx','vh','vi','v')

subplot(3,4,3); hold on; box on
plot(Topt,T, MarkerSize=3)
plot(Topt,T.*sin(theta), MarkerSize=3)
plot(Topt,T.*cos(theta), MarkerSize=3)
xlabel('t'); ylabel('T');
legend('T','Tx','Th',Location='south')

subplot(3,4,4); hold on; box on
plot(Topt,rad2deg(theta),LineWidth=2)
plot(Topt,rad2deg(gamma))
plot(Topt,rad2deg(alpha))
xlabel('t'); ylabel('theta');
legend('theta','gamma','alpha',Location='south')

subplot(3,4,5); hold on; box on
plot(Topt,gradient(vx,Topt), '-o', MarkerSize=3)
plot(Topt,gradient(vh,Topt), '-o', MarkerSize=3)
xlabel('t'); ylabel('a');
legend('ax','ah',Location='south')

subplot(3,4,6); hold on; box on
plot(Topt,Pinduced)
plot(Topt,Pclimb)
plot(Topt,Pprofile)
plot(Topt,Power,'LineWidth',2)
Popt = [Pinduced,Pclimb,Pprofile,Power];
xlabel('t'); ylabel('P');

subplot(3,4,7); hold on; box on
plot(Topt,SOC, '-o', MarkerSize=3)
xlabel('t'); ylabel('SOC');

subplot(3,4,8); hold on; box on
plot(Topt,Up, '-o', MarkerSize=3)
plot(Topt,Uoc, '-o', MarkerSize=3)
plot(Topt,U, '-o', MarkerSize=3)
xlabel('t'); ylabel('U');
legend('Up','Uoc','U',Location='south')

subplot(3,4,9); hold on; box on
plot(Topt,I)
plot(Topt,Xopt(:,8))
xlabel('t'); ylabel('I');

subplot(3,4,10); hold on; box on
plot(Topt,Ro, '-o', MarkerSize=3)
plot(Topt,Rp, '-o', MarkerSize=3)
xlabel('t'); ylabel('R');

subplot(3,4,11); hold on; box on
plot(Topt,Cp, '-o', MarkerSize=3)
xlabel('t'); ylabel('Cp');

subplot(3,4,12); hold on; box on
plot(Topt,TempC, '-o', MarkerSize=3)
xlabel('t'); ylabel('Temp');
%%
% energy = sum(Power)/1000/3600;
% fprintf('Energy = %s /kWh\n',num2str(energy*params.nBat));
save('NPplot.mat','Power','Crate','TempC','Topt','SOC')
end

