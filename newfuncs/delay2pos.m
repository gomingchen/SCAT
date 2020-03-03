
function d = delay2pos(coords_ears, delays)
% calculate the coordinates of ONE glints with four delays in ms from two
% ears
%%% ATTENTION! the delays are doubled in the recording

% left ear, center of left circle
x1 = coords_ears(1,1);
y1 = coords_ears(1,2);

% right ear, center of right circle
x2 = coords_ears(2,1);
y2 = coords_ears(2,2);

% center of the two origins

R = dista(coords_ears(1,:), coords_ears(2,:)); % distance between two centers of the circles

xc = (x1 + x2)/2;
yc = (y1 + y2)/2;

% the radii of left circle
r = delays*1E-6/2*340; % radius in meters
r1 = r(1);
r2 = r(2);

a = r1^2 + r2^2;
b = r1^2 - r2^2;

xg1 = xc + b/(2*R^2)*(x2 - x1) + 0.5*sqrt(2*a/R^2 - (b/R^2)^2 - 1) * (y2 - y1);
xg2 = xc + b/(2*R^2)*(x2 - x1) - 0.5*sqrt(2*a/R^2 - (b/R^2)^2 - 1) * (y2 - y1);

yg1 = yc + b/(2*R^2)*(y2 - y1) + 0.5*sqrt(2*a/R^2 - (b/R^2)^2 - 1) * (x1 - x2);
yg2 = yc + b/(2*R^2)*(y2 - y1) - 0.5*sqrt(2*a/R^2 - (b/R^2)^2 - 1) * (x1 - x2);

d = [xg1 yg1 xg2 yg2];
end

function dd = dista(a,b)
    dd = sqrt((a(1) - b(1))^2 + (a(2) - b(2))^2);
end