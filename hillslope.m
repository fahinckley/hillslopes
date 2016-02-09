%------------------------------
% Hillslope modeling
%------------------------------
% Franklin Hinckley
% 8 February 2016
%------------------------------

%% Clean up workspace
clearvars
close all
clc

%% Initialize hillslope
dx = 0.25; % position step [m]
x = -50:dx:50; % position vector [m]
z0 = 10; % intial surface height [m]
z = z0*ones(size(x)); % uniform surface
zR = 0; % height of rock interface

%% Define edge condition and efficiency
eDot = -0.02; % incision rate [m/yr]
k = 250; % transport efficiency [kg/yr]
zstar = 3; % characteristic depth [m]
rhoR = 2400; % density of rock [kg/m^3]
rhoS = 1500; % density of soil [kh/m^3]
wdot = -0.01; % weathering rate [m/yr]

%% Set up time vector
dt = 0.1; % time step [yr]
t = 0:dt:500;

%% Main loop
zS = zeros(length(z),length(t));
zS(:,1) = z;
for ii = 1:length(t)
    
    % Compute slope
    dzdx = diff(z)/dx;

    % Get box-centered heights
    zM = zeros(size(dzdx)); % box-centered surface height
    for jj = 1:length(dzdx)
        zM(jj) = mean(z(jj:jj+1));
    end
    
    % Compute flux
    q = -k.*dzdx;%.*(1-exp(-(zM-zR)/zstar));

    % Compute height rate
    dqdx = diff(q)/dx; 
    dHdt = (rhoR/rhoS)*wdot - (1/rhoS)*dqdx;
    
    % Add edge condition
    dHdt = [eDot dHdt eDot];

    % Update height
    z = z + dHdt*dt;

    % Update rock boundary
    zR = zR + wdot*dt;
    
    % Save output
    zS(:,ii) = z;
end

%% Plots
% figure
% hold on
% plot(x,zS')
% plot(boundLow,z,'-k','Linewidth',2)
% plot(boundHigh,z,'-k','Linewidth',2)
% plot([0 0],[min(z) max(z)],'--r')
% plot([min(T(:))-2 max(T(:))+2],[z(zInd) z(zInd)],'--k')
% hold off
% xlim([min(boundLow)-2 max(boundHigh)+2])
% ylim([0 zmax])
% set(gca,'Ydir','reverse')
% title('Geotherm','Fontsize',14)
% ylabel('Depth [m]','Fontsize',12)
% xlabel('Temp [\circC]','Fontsize',12)

% Animation
figure
for ii = 1:length(t)
    plot(x,zS(:,ii),'b')
    hold on
%     plot(boundLow,z,'-k')
%     plot(boundHigh,z,'-k')
%     plot([0 0],[min(z) max(z)],'r')
    hold off
%     xlim([min(TS(:))-2 max(TS(:))+2])
%     ylim([0 zmax])
%     set(gca,'Ydir','reverse')
%     title('Geotherm','Fontsize',14)
%     ylabel('Depth [m]','Fontsize',12)
%     xlabel('Temp [\circC]','Fontsize',12)
%    M(:,ii) = getframe(gcf);
    pause(0.01)    
end
