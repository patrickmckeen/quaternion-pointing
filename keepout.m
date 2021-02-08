%Example calculation of possible eigenaxis rotational trajectories to align
%two vectors.
%includes shortest slew, as well as avoiding a "keep-out" zone
%Copyright Patrick McKeen, February 2021

%% INPUTS
u = [0 1 -1].'/sqrt(2); %desired target vector to point at in rotated world frame
v = [1 0 0].'; %desired body vector to point in body frame
s = [1 sqrt(3)/2 1/2].'/sqrt(2); %vector towards sun in rotated world frame
t = [0 1 0].'; %imager boresight vector in body frame
N = 1000; %values for alpha (between 0 and 2*pi) to consider
eta = pi/6; %half-cone angle of "keep-out" zone in radians
minscale = 100;

%% CALCULATION
u = u/norm(u);
v = v/norm(v);
s = s/norm(s);
t = t/norm(t);
alpha = linspace(0,2*pi,N);
theta = acos(dot(v,u));
x = cross(v,u)/norm(cross(v,u));
y = (v+u)/norm(v+u);
k = x*cos(alpha) + y*sin(alpha);
gamma = pi-2*atan(sin(alpha)*cot(theta/2));
o = ones(1,N);
cosEta = cos(eta);
alphaMin = 0; %by definition if q = [1 0 0 0]

stDot = dot(s*o,t*o);
stCrossKDot = dot(s*o,cross(t*o,k));
sCrossKTCrossKDot = dot(cross(s*o,k),cross(t*o,k));
m = sqrt(sCrossKTCrossKDot.^2 + stCrossKDot.^2);

rhoZeroLim = stDot - cosEta;
rhoGammaLim = stDot - 2*((sCrossKTCrossKDot + stCrossKDot.*sin(alpha)*cot(theta/2))./(1 + sin(alpha).^2*cot(theta/2)^2)) - cosEta;
rhoMaxVal0 = gamma - mod(atan2(stCrossKDot,sCrossKTCrossKDot),2*pi);
rhoMaxValTest = tan(theta/2)*((m + sCrossKTCrossKDot)./(stCrossKDot)) - sin(alpha);
rhoMaxLim = stDot - sCrossKTCrossKDot + m - cosEta;
viableAlpha = all([rhoZeroLim < 0; rhoGammaLim < 0; or((rhoMaxLim < 0),(rhoMaxValTest < 0))],1);

%% PLOT
figure;
e0 = plot(alpha,0*alpha,'k');
hold on
e1 = plot(alpha,rhoZeroLim/min(max(abs(rhoZeroLim)),minscale),'b');
e2 = plot(alpha,rhoGammaLim/min(max(abs(rhoGammaLim)),minscale),'r');
e3 = plot(alpha,rhoMaxLim/min(max(abs(rhoMaxLim)),minscale),'-.');
e4 = plot(alpha,rhoMaxValTest/min(max(abs(rhoMaxValTest)),minscale),'m-.');
e5 = plot(alpha(viableAlpha),0*alpha(viableAlpha),'go');
e6 = plot(alphaMin,0,'*');
axis([0 2*pi, -1.5,1.5])
xlabel('\alpha')
ylabel('Limit Values')
legend([e1,e2,e3,e4,e5,e6],'\rho=0','\rho=\gamma','\rho=\rho_{max}','\rho_{max}<\gamma','Acceptable \alpha','\alpha with minimum slew')

%% RESULTS
acceptableAlpha = alpha(viableAlpha);
alphaRising = alpha([false and(~viableAlpha(1:end-1),viableAlpha(2:end))]);
alphaFalling = alpha([and(viableAlpha(1:end-1),~viableAlpha(2:end)) false]);

if viableAlpha(1)
    alphaRising = [0 alphaRising];
end


if viableAlpha(end)
    alphaFalling = [alphaFalling 2*pi];
end

acceptableAlphaRanges = [alphaRising.' alphaFalling.'];

sprintf('Acceptable alphas are in the following ranges:')
for i = 1:length(alphaRising)
    sprintf('%.3f <= alpha <= %.3f',alphaRising(i),alphaFalling(i))
end