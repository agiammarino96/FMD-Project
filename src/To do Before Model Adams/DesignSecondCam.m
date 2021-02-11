%% Preliminary Study
close all
clear all
clc

% Initial and Final Pressure Angle as beta_0 and link length change
M=80; % Variation beta_0
U=10; % Variation link length
press_ang_initial = zeros(U,M);
press_ang_final = zeros(U,M);
beta_0_var = zeros(1,M);
link_var = zeros(1,U);
beta_0_min = zeros(1,U);
index_beta_0_min = zeros(1,U);
for k=1:U
    for j=1:M
        link=6 + k*2; % length of the link
        link_var(k) = link;
        Rr=10; % roller radius
        beta_0=1*j*pi/180; % initial slope of the link with the vertical axis
        beta_0_var(j) = beta_0*180/pi;
        
        load('data_x.mat'); % this is the motiolaw of x obtained from TotalDesignFirstCamSymmetric
        
        load('data_theta.mat'); % this is the motiolaw of theta obtained from motionlaw_Theta
        
        beta=beta_0+data_theta.pos;
        v_t=link*data_theta.vel;
        v_a=sqrt((v_t.*cos(beta)+data_x.vel).^2+(v_t.*sin(beta)).^2);
        
        gamma=acos((v_t.^2+v_a.^2-data_x.vel.^2)./(2*v_t.*v_a));
        n=find(data_theta.pos~=0,1); % this is the point where beta starts to move
        gamma(1:n-1)=beta_0; % at the start v_t=0 because theta_vel=0 and would give NaN but in this points theta=beta_0
        gamma(end)=gamma(end-1);
        
        press_ang=pi/2-gamma;
        press_ang_initial(k,j) = press_ang(1)*180/pi;
        press_ang_final(k,j) = press_ang(end)*180/pi;
        
        a_n=data_theta.vel.^2*link.*cos(gamma)+data_theta.acc*link.*sin(gamma)-data_x.acc.*cos(pi/2-beta+gamma);
        rho0=(v_a.^2)./a_n;
        
        r_i=sqrt(link^2+Rr^2-2*link*Rr*cos(gamma));
        phi_i=beta+asin((Rr*sin(gamma))./r_i);
        
        r_e=sqrt(link^2+Rr^2+2*link*Rr*cos(gamma));
        phi_e=beta-asin((Rr*sin(gamma))./r_e);
        
        
        [x_pitch,y_pitch]=pol2cart(3*pi/2+beta,link);
        x_pitch=-(x_pitch+data_x.pos);
        
        [x_i,y_i]=pol2cart(3*pi/2+phi_i,r_i);
        x_i=-(x_i+data_x.pos);
        
        [x_e,y_e]=pol2cart(3*pi/2+phi_e,r_e);
        x_e=-(x_e+data_x.pos);
    end
    [beta_0_min(k), index_beta_0_min(k)] = min(press_ang_final(k,:));
end

% Plots pressure angle study

figure
for k=1:U
    plot(beta_0_var,press_ang_final(k,:),'LineWidth',2,'DisplayName',['link length [mm] = ', num2str(link_var(k))])
    hold on
    plot(beta_0_var(index_beta_0_min(k)),beta_0_min(k),'o','HandleVisibility','off')
    grid on
end
title('Final Pressure Angle')
%legend(['link length [mm] = ', num2str(link_var(1))],['link length [mm] = ', num2str(link_var(2))],['link length [mm] = ', num2str(link_var(3))],['link length [mm] = ', num2str(link_var(4))],['link length [mm] = ', num2str(link_var(5))])
legend('show')
xlabel('\beta_0 [°]')
ylabel('\theta [°]')

figure
for k=1:U
    plot(beta_0_var,press_ang_initial(k,:),'LineWidth',2,'DisplayName',['link length [mm] = ', num2str(link_var(k))])
    hold on
    grid on
end
title('Initial Pressure Angle')
%legend(['link length [mm] = ', num2str(link_var(1))],['link length [mm] = ', num2str(link_var(2))],['link length [mm] = ', num2str(link_var(3))],['link length [mm] = ', num2str(link_var(4))],['link length [mm] = ', num2str(link_var(5))])
legend('show')
xlabel('\beta_0 [°]')
ylabel('\theta [°]')
%% Generation cam


% First of all run TotalDesignFirstCamSymmetric to have the motionlaw of x
% and motionlaw_Theta to have the motionlaw of theta


link=20; % length of the link
Rr=10; % roller radius
beta_0=23*pi/180; % initial slope of the link with the vertical axis

load('data_x.mat'); % this is the motiolaw of x obtained from TotalDesignFirstCamSymmetric

load('data_theta.mat'); % this is the motiolaw of theta obtained from motionlaw_Theta

beta=beta_0+data_theta.pos;
v_t=link*data_theta.vel;
v_a=sqrt((v_t.*cos(beta)+data_x.vel).^2+(v_t.*sin(beta)).^2);

gamma=acos((v_t.^2+v_a.^2-data_x.vel.^2)./(2*v_t.*v_a));
n=find(data_theta.pos~=0,1); % this is the point where beta starts to move
gamma(1:n-1)=beta_0; % at the start v_t=0 because theta_vel=0 and would give NaN but in this points theta=beta_0
gamma(end)=gamma(end-1);

press_ang=pi/2-gamma;

a_n=data_theta.vel.^2*link.*cos(gamma)+data_theta.acc*link.*sin(gamma)-data_x.acc.*cos(pi/2-beta+gamma);
rho0=(v_a.^2)./a_n;

r_i=sqrt(link^2+Rr^2-2*link*Rr*cos(gamma));
phi_i=beta+asin((Rr*sin(gamma))./r_i);

r_e=sqrt(link^2+Rr^2+2*link*Rr*cos(gamma));
phi_e=beta-asin((Rr*sin(gamma))./r_e);


[x_pitch,y_pitch]=pol2cart(3*pi/2+beta,link);
x_pitch=-(x_pitch+data_x.pos);

[x_i,y_i]=pol2cart(3*pi/2+phi_i,r_i);
x_i=-(x_i+data_x.pos);

[x_e,y_e]=pol2cart(3*pi/2+phi_e,r_e);
x_e=-(x_e+data_x.pos);


%% close the profile
N=length(data_x.pos);
eps=3;

xm_right=x_pitch(1)+eps;
ym_right=y_pitch(1);
alpha_right=linspace(-pi/2,pi/2,N);
semicircle_right_x=xm_right+Rr*cos(alpha_right);
semicircle_right_y=ym_right+Rr*sin(alpha_right);

angular_coeff=(y_i(end)-y_e(end))/(x_i(end)-x_e(end));
angle_left=atan(angular_coeff);
xm_left=x_pitch(end)+eps*cos(angle_left+pi/2);
ym_left=y_pitch(end)+eps*sin(angle_left+pi/2);
alpha_left=linspace(angle_left,angle_left+pi,N);
semicircle_left_x=xm_left+Rr*cos(alpha_left);
semicircle_left_y=ym_left+Rr*sin(alpha_left);



profileCam_x(1,1:N) = semicircle_left_x;
profileCam_x(1,(N+1):2*N) = flip(x_e);
profileCam_x(1,(2*N+1):3*N) = semicircle_right_x;
profileCam_x(1,(3*N+1):4*N) = x_i;
profileCam_x(1)=profileCam_x(end);


profileCam_y(1,1:N) = semicircle_left_y;
profileCam_y(1,(N+1):2*N) = flip(y_e);
profileCam_y(1,(2*N+1):3*N) = semicircle_right_y;
profileCam_y(1,(3*N+1):4*N) = y_i;
profileCam_y(1)=profileCam_y(end);



%% plots of cam profile

figure
plot(profileCam_x,profileCam_y,'Displayname','Cam Profile')
hold on
% plot(profileCam_x,profileCam_y,'o','Displayname','Cam Profile')
plot(x_pitch,y_pitch,'Displayname','Pitch Profile')
grid on
legend
xlabel('X [mm]')
ylabel('Y [mm]')
title('Profile')
xlim ([-300 10]) %these 3 last rows of code are to center the plot and keep the same axis proportion
ylim ([-150 160])
daspect([1 1 1])

% figure
% plot(x_pitch,y_pitch,'Displayname','Pitch Profile')
% hold on
% plot(x_i,y_i,'Displayname','Internal Cam Profile')
% hold on
% plot(x_e,y_e,'Displayname','External Cam Profile')
% grid on
% legend
% xlabel('X [mm]')
% ylabel('Y [mm]')
% title('Profile')
% xlim ([-300 10]) %these 3 last rows of code are to center the plot and keep the same axis proportion
% ylim ([-150 160])
% daspect([1 1 1])

%% exporting the profile

profile = [profileCam_x ; profileCam_y; zeros(1,4*N)]';
pitch = [x_pitch ; y_pitch ; zeros(1,length(x_pitch))]';
profileTab = table(profile);
pitchTab = table(pitch);
writetable(profileTab,'SecondCam.txt');
writetable(pitchTab,'pitchProfileSecondCam.txt')


%% pressure angle

Abscissa = [100 30]; % we design the profile using only the rise
abscissa = linspace (0 ,sum( Abscissa ) ,N); % vector with N elements

figure
plot(abscissa,(press_ang)*180'/pi,'Displayname','Pressure Angle') %in deg
title('Pressure Angle')
xlabel('\alpha [°]')
ylabel('\theta [°]')

%% undercut check
k=0;
for i=1:N
    if(abs(rho0(i))<Rr)
        k=k+1;
    end
end
k

m=find(abs(rho0)<Rr) % these are the points where we have undercut, if we have it only in the last point I think it is negligible

%% Curvature radius for contact pressure computation

rho = abs(rho0)-Rr;
curvature_radii_secondCam.rho = rho;
curvature_radii_secondCam.rho0 = rho0;

save('curvature_radii_secondCam.mat','curvature_radii_secondCam');

%% Stuff for contact force computation second cam along X
% %% Trend of beta
% figure
% plot(abscissa,beta*180/pi)
% title('Angle Follower 2')
% xlabel('\alpha [°]')
% ylabel('\beta [°]')
% 
% %% Study of Forces for Spring design
% % Here we assume X positive leftwards and negative rightwards
% % With this convention the force which detouches the first follower from
% % the cam is the negative one
% 
% % Initial data needed
% 
% % data_x contains the rise of the first motion law (geometrical)
% % data_theta contains the rise of the second motion law (geometrical)
% % data_x.vel [mm/rad], data_theta.vel [rad/rad]
% m_tot = 4.5 + 1.2; % [kg]
% J_tot = 0.5; % [kg*m^2]
% omega = pi; % [rad/s]
% h = 180; % [mm]
% h_meters = h*10^(-3); % [m]
% link_meters = link*10^(-3); % [m]
% acc_time_x = data_x.acc*omega^2*10^(-3); % [m/s^2]
% acc_time_theta = data_theta.acc*omega^2; % [rad/s^2]
% pos_time_x = data_x.pos*10^(-3); % [m]
% pos_time_theta = data_theta.pos; % [rad]
% t_final_rise = 130*pi/(180*omega);
% t_rise = linspace(0,t_final_rise,N); % [s]
% 
% % Important plots (position and acceleration) motion laws
% figure
% yyaxis right
% plot(t_rise,pos_time_x)
% ylabel('position [m]')
% grid on
% yyaxis left
% plot(t_rise,acc_time_x)
% ylabel('acceleration [m/s^2]')
% title('Position and Acceleration First Follower Rise')
% xlabel('time [s]')
% xlim([t_rise(1) t_rise(end)])
% 
% figure
% yyaxis left
% plot(t_rise,pos_time_theta)
% ylabel('position [rad]')
% grid on
% yyaxis right
% plot(t_rise,acc_time_theta)
% ylabel('acceleration [rad/s^2]')
% title('Position and Acceleration Second Follower Rise')
% xlabel('time [s]')
% xlim([t_rise(1) t_rise(end)])
% 
% 
% % Computation Inertial forces along X: due to translating mass
% Fi_x = -acc_time_x*m_tot; % [N]
% 
% % Computation Inertial forces along X: due to rotating follower
% Fi_theta = -acc_time_theta*J_tot/link_meters; % [N]
% S_theta = Fi_theta./cos(press_ang); 
% zeta = press_ang + beta - pi/2;
% S_theta_x = S_theta.*sin(zeta);
% S_theta_y = S_theta.*cos(zeta);
% 
% % Plots
% figure
% plot(t_rise,Fi_x)
% grid on
% xlabel('time [s]')
% ylabel('F_i translation [N]')
% xlim([t_rise(1) t_rise(end)])
% title('Load due to translation')
% 
% figure
% plot(t_rise,S_theta_x)
% grid on
% xlabel('time [s]')
% ylabel('F_i rotation [N]')
% xlim([t_rise(1) t_rise(end)])
% title('Load due to rotation')
% 
% figure
% plot(t_rise,press_ang)
% hold on
% grid on
% plot(t_rise,beta)
% plot(t_rise,zeta)
% legend('Pressure Angle','\beta','\zeta')
% xlabel('time [s]')
% ylabel('Angle [rad]')
% xlim([t_rise(1) t_rise(end)])
% 
% 
% % Spring design
% 
% totalLoad = S_theta_x+Fi_x;
% k = 0; %[N/m]
% preload = -min(totalLoad); %[N]
% forceDisplacement = +pos_time_x*k-preload;
% 
% figure
% plot(t_rise,forceDisplacement)
% hold on
% plot(t_rise,totalLoad)




