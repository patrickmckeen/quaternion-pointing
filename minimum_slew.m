%Example calculation of minimum slew
%Copyright Patrick McKeen, February 2021

%% INPUTS
u = [0 1 -1].'/sqrt(2); %desired target vector to point at in rotated world frame
v = [1 0 0].'; %desired body vector to point in body frame
p = [4 1 -2 2].'/5;%starting orientation quaternion

%% CALCULATION
u = u/norm(u);
v = v/norm(v);
theta = acos(dot(v,u));
x = cross(v,u)/norm(cross(v,u));
y = (v + u)/norm(v + u);
xquat = [cos(theta/2); x*sin(theta/2)];
yquat = [0;y];
beta1 = atan(xquat.'*p/(yquat.'*p));
beta2 = beta1 + pi;
f1 = ((xquat.'*p)*xquat + (yquat.'*p)*yquat)/sqrt((xquat.'*p)^2 + (yquat.'*p)^2);
f2 = -f1;

pinv = [p(1);-p(2:end)];

rot1 = [pinv(1)*f1(1) - dot(pinv(2:end),f1(2:end));...
        f1(1)*pinv(2:end) + pinv(1)*f1(2:end) + cross(pinv(2:end),f1(2:end))];
rot2 = [pinv(1)*f2(1) - dot(pinv(2:end),f2(2:end));...
        f2(1)*pinv(2:end) + pinv(1)*f2(2:end) + cross(pinv(2:end),f2(2:end))];

ax1 = rot1(2:end)./norm( rot1(2:end));
ang1 = 2*atan2(norm( rot1(2:end)),rot1(1));
ax2 = rot2(2:end)./norm( rot2(2:end));
ang2 = 2*atan2(norm(rot2(2:end)),rot2(1));
    
%% RESULTS

sprintf('Minimum slew rotation:')
sprintf('Rotation')
sprintf('-- Goal Orientation Quaternion: [ %.3f %.3f %.3f %.3f ]^T',f1)
sprintf('-- Corresponding Beta: %.3f radians',beta1)
sprintf('-- Rotation from Starting Orientation (Quaternion): [ %.3f %.3f %.3f %.3f ]^T',rot1)
sprintf('-- Axis of Rotation (Body Frame): [ %.3f %.3f %.3f ]^T ',ax1)
sprintf('-- Angle of Rotation (Body Frame): %.3f radians',ang1)

%% NOT RELEVANT--shows that using f2=-f1 gives the same results, 
% sprintf('Rotation 2')
% sprintf(['-- Goal Orientation Quaternion:' num2str(f2.')])
% sprintf(['-- Corresponding Beta:' num2str(beta2.')])
% sprintf(['-- Rotation from Starting Orientation:' num2str(rot2.')])
% sprintf(['-- Axix of Rotation (Body Frame):' num2str(ax2.')])
% sprintf(['-- Angle of Rotation (Body Frame):' num2str(ang2.')])