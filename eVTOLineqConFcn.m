function cineq = eVTOLineqConFcn(stage,x,u,dmv,p)

persistent obstacles eVTOL
eWidth = p(1+68);
eHeight = p(2+68);
ObsX = p(5+68);
ObsH = p(6+68);
ObsWidth = p(7+68);
ObsHeight = p(8+68);

safetyDistance = 10;
if isempty(obstacles)
    %% eVTOL
    eVTOL = mpcShape('Rectangle', eWidth, eHeight);
    %% Obstacles
    obs1 = mpcShape('Rectangle',ObsWidth,ObsHeight);
    [obs1.X, obs1.Y] = deal(ObsX,ObsH);
    obstacles = {obs1};
end

% Update eVTOL
eVTOL.X = x(3);
eVTOL.Y = x(4);
eVTOL.Theta = u(2);

% Calculate distances from trailer to obstacles
numObstacles = numel(obstacles);
distance = zeros(numObstacles,1);
for ct = 1:numObstacles
    [collisionStatus1, distance(ct)] = ...
        controllib.internal.gjk.Base2d.checkCollision(eVTOL, obstacles{ct});

    if collisionStatus1
        distance(ct) = -10^4;
    end
end
cineq = -distance + safetyDistance;