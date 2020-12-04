%% RISE MOTION LAW FIRST CAM %%%%%%%%%%%%
close all
clear all
clc

% Choose d1, d2, d3 and d4; compute d5, d6 and d7
d = zeros(1,7);
d(1) = 0.2;
d(2) = 0.1;
d(3) = 0.25;
d(4) = 0.0;
d(5) = 1-(30/130+sum(d(1:4))); % chosen so that the sum(d(1:5)) is so that alpha_5=100° (d6+d7=30/130)
h = 180; % Total rise motion law
N = 100; % Number of samples
Ab = 130; % Cam angle correspondent to rise
t_geom1 = linspace (0, Ab, N); % Cam angle discretization
t_adim1 = linspace(0,1,N); % Discretization for adimensional motion law

% Computation of A, B, d2 and d6 so that
% Parametrization of boundary conditions for each segment
syms A B d6 d7
y11 = A*2*d(1)/pi;  %y11 stands for y1'
y1 =  A*2*d(1)^2/pi*(1-2/pi);
y21 = y11 +A*d(2);
y2 = y1+A*d(2)*(2*d(1)/pi+d(2)/2);
y31 = y21 +A*2*d(3)/pi;
y3 = y2+A*d(3)*(2*d(1)/pi+d(2)+4*d(3)/(pi^2));
y41 = y31;
y4 = y3+A*d(4)*(2*d(1)/pi+d(2)+2*d(3)/pi);
y51 = y41-B*2*d(5)/pi;
y5 = y4+A*d(5)*(2*d(1)/pi+d(2)+2*d(3)/pi)-B*2*d(5)^2/pi*(1-2/pi);
y61 = y51-B*d6;
y6 = y5+y51*d6-B*d6^2/2;
y71 = y61-B*2*d7/pi;
y7 = y6+y61*d7-B*(2*d7/pi)^2;

% Conditions to be satisfied:
% 1. y7==1 means that at the end of the motion law we want y=180mm
% 2. y71==0 means that at the end of the motion law we want zero velocity
% 3. y5==160/180 means that at point 5 of the motion law we want y=160mm
% 4. d6+d7==30/130 means that at point 5 we want alpha=100°

[A, B, d6, d7] = solve(y7==1,y71==0,y5==160/180,d6+d7==30/130);

% Values of A, B and d computed converted to double
A_solution = double(A);
B_solution = double(B);
d6_solution = double(d6);
d7_solution = double(d7);

% Choice of the solution to use
A = A_solution(2);
B = B_solution(2);
d(6) = d6_solution(2);
d(7) = d7_solution(2);

% Computation of motion law
y11 = A*2*d(1)/pi;  %y11 stands for y1'
y1 =  A*2*d(1)^2/pi*(1-2/pi);
y21 = y11 +A*d(2);
y2 = y1+A*d(2)*(2*d(1)/pi+d(2)/2);
y31 = y21 +A*2*d(3)/pi;
y3 = y2+A*d(3)*(2*d(1)/pi+d(2)+4*d(3)/(pi^2));
y41 = y31;
y4 = y3+A*d(4)*(2*d(1)/pi+d(2)+2*d(3)/pi);
y51 = y41-B*2*d(5)/pi;
y5 = y4+A*d(5)*(2*d(1)/pi+d(2)+2*d(3)/pi)-B*2*d(5)^2/pi*(1-2/pi);
y61 = y51-B*d(6);
y6 = y5+y51*d(6)-B*d(6)^2/2;
y71 = y61-B*2*d(7)/pi;
y7 = y6+y61*d(7)-B*(2*d(7)/pi)^2;

pos_adim1 = zeros (1,N);
vel_adim1 = zeros (1,N);
acc_adim1 = zeros (1,N);
jerk_adim1 = zeros(1,N);
pos_geom1 = zeros (1,N);
vel_geom1 = zeros (1,N);
acc_geom1 = zeros (1,N);
jerk_geom1 = zeros(1,N);

for i=1:length(d)
    a(i)=sum(d(1:i));
end

for k=1:N
    epsilon = t_geom1(k)/Ab;
    
    if 0 <= epsilon && epsilon <=d(1)
        out.j = A*pi/(2*d(1))*cos(epsilon*pi/(2*d(1)));
        out.a =  A*sin(epsilon*pi/(2*d(1)));
        out.v =  A*2/pi*d(1)*(1-cos(epsilon*pi/(2*d(1))));
        out.s =  A*2/pi*d(1)*(epsilon-2*d(1)/pi*sin(epsilon*pi/(2*d(1))));
        
        
    elseif d(1) < epsilon && epsilon <=sum(d(1:2))
        out.j = 0;
        out.a = A;
        out.v = y11+A*(epsilon-a(1));
        out.s = y1+y11*(epsilon-a(1))+A*(epsilon-a(1))^2/2;
        
        
    elseif sum(d(1:2)) < epsilon && epsilon <=sum(d(1:3))
        out.j = -A*pi/(2*d(3))*sin((epsilon-a(2))*pi/(2*d(3)));
        out.a = A*cos((epsilon-a(2))*pi/(2*d(3)));
        out.v = y21+A*2*d(3)/pi*sin((epsilon-a(2))*pi/(2*d(3)));
        out.s = y2+y21*(epsilon-a(2))+A*(2*d(3)/pi)^2*(1-cos((epsilon-a(2))*pi/(2*d(3))));
        
        
    elseif sum(d(1:3)) < epsilon && epsilon <=sum(d(1:4))
        out.j = 0;
        out.a = 0;
        out.v = y31;
        out.s = y3+y31*(epsilon-a(3));
        
        
    elseif sum(d(1:4)) < epsilon && epsilon <=sum(d(1:5))
        out.j = -B*pi/(2*d(5))*cos((epsilon-a(4))*pi/(2*d(5)));
        out.a = -B*sin((epsilon-a(4))*pi/(2*d(5)));
        out.v = y41-B*2*d(5)/pi*(1-cos((epsilon-a(4))*pi/(2*d(5))));
        out.s = y4+y41*(epsilon-a(4))-B*2*d(5)/pi...
            *(epsilon-a(4)-2*d(5)/pi*sin((epsilon-a(4))*pi/(2*d(5))));
        
        
    elseif sum(d(1:5)) < epsilon && epsilon <=sum(d(1:6))
        out.j = 0;
        out.a = -B;
        out.v = y51-B*(epsilon-a(5));
        out.s = y5+y51*(epsilon-a(5))-B*(epsilon-a(5))^2/2;
        
        
    else
        out.j = B*pi/(2*d(7))*sin((epsilon-a(6))*pi/(2*d(7)));
        out.a = -B*cos((epsilon-a(6))*pi/(2*d(7)));
        out.v = y61-B*2*d(7)/pi*sin((epsilon-a(6))*pi/(2*d(7)));
        out.s = y6+y61*(epsilon-a(6))-B*(2*d(7)/pi)^2*(1-cos((epsilon-a(6))*pi/(2*d(7))));
        
    end
    pos_geom1(k) = out.s*h;
    vel_geom1(k) = out.v*h/Ab;
    acc_geom1(k) = out.a*h/Ab^2;
    jerk_geom1(k) = out.j*h/Ab^3;
    pos_adim1(k) = out.s;
    vel_adim1(k) = out.v;
    acc_adim1(k) = out.a;
    jerk_adim1(k) = out.j;
end


%% DWELL MOTION LAW FIRST CAM %%%%%%%%%%%%%%%%%%

% Dwell between 130° and 170°
t_geom_dwell1 = linspace(0,170-130,N);
pos_geom_dwell1 = zeros(1,N);
vel_geom_dwell1 = zeros(1,N);
acc_geom_dwell1 = zeros(1,N);
jerk_geom_dwell1 = zeros(1,N);

% Dwell between 300° and 360°
t_geom_dwell2 = linspace(0,360-300,N);
pos_geom_dwell2 = zeros(1,N);
vel_geom_dwell2 = zeros(1,N);
acc_geom_dwell2 = zeros(1,N);
jerk_geom_dwell2 = zeros(1,N);

% Dwell between 130° and 170°
t_adim_dwell1 = linspace(0,(170-130)/360,N);
pos_adim_dwell1 = zeros(1,N);
vel_adim_dwell1 = zeros(1,N);
acc_adim_dwell1 = zeros(1,N);
jerk_adim_dwell1 = zeros(1,N);

% Dwell between 300° and 360°
t_adim_dwell2 = linspace(0,(360-300)/360,N);
pos_adim_dwell2 = zeros(1,N);
vel_adim_dwell2 = zeros(1,N);
acc_adim_dwell2 = zeros(1,N);
jerk_adim_dwell2 = zeros(1,N);


%% FALL MOTION LAW FIRST CAM %%%%%%%%%%%%
t_geom2 = linspace (0, Ab, N); % Cam angle discretized
t_adim2 = linspace(0,1,N); % Discretization for adimensional law


pos_adim2 = zeros (1,N);
vel_adim2 = zeros (1,N);
acc_adim2 = zeros (1,N);
jerk_adim2 = zeros(1,N);
pos_geom2 = zeros (1,N);
vel_geom2 = zeros (1,N);
acc_geom2 = zeros (1,N);
jerk_geom2 = zeros(1,N);

for j=0:N-1
    pos_geom2(j+1) = pos_geom1(end-j);
    vel_geom2(j+1) = -vel_geom1(end-j);
    acc_geom2(j+1) = acc_geom1(end-j);
    jerk_geom2(j+1) = -jerk_geom1(end-j);
    pos_adim2(j+1) = pos_adim1(end-j);
    vel_adim2(j+1) = -vel_adim1(end-j);
    acc_adim2(j+1) = acc_adim1(end-j);
    jerk_adim2(j+1) = -jerk_adim1(end-j);
end


%% TOTAL LAW %%%%%%%%%%%%%%%%%%%%%%%%
% ADIMENSIONAL
% Allocation of space
t_adim_tot = zeros(1,4*N);
pos_adim_tot = zeros(1,4*N);
vel_adim_tot = zeros(1,4*N);
acc_adim_tot = zeros(1,4*N);
jerk_adim_tot = zeros(1,4*N);

% alpha
t_adim_tot(1:N) = linspace(0,130/360,N);
t_adim_tot(N+1:2*N) = linspace(130/360,170/360,N);
t_adim_tot(2*N+1:3*N) = linspace(170/360,300/360,N);
t_adim_tot(3*N+1:4*N) = linspace(300/360,1,N);

% position
pos_adim_tot(1:N) = pos_adim1;
pos_adim_tot(N+1:2*N) = pos_adim1(end)+pos_adim_dwell1;
pos_adim_tot(2*N+1:3*N) = pos_adim2;
pos_adim_tot(3*N+1:4*N) = pos_adim_tot(3*N)+pos_adim_dwell2;

% velocity
vel_adim_tot(1:N) = vel_adim1;
vel_adim_tot(N+1:2*N) = vel_adim1(end)+vel_adim_dwell1;
vel_adim_tot(2*N+1:3*N) = vel_adim2;
vel_adim_tot(3*N+1:4*N) = vel_adim_tot(3*N)+vel_adim_dwell2;

% acceleration
acc_adim_tot(1:N) = acc_adim1;
acc_adim_tot(N+1:2*N) = acc_adim1(end)+acc_adim_dwell1;
acc_adim_tot(2*N+1:3*N) = acc_adim2;
acc_adim_tot(3*N+1:4*N) = acc_adim_tot(3*N)+acc_adim_dwell2;

% jerk
jerk_adim_tot(1:N) = jerk_adim1;
jerk_adim_tot(N+1:2*N) = jerk_adim_dwell1;
jerk_adim_tot(2*N+1:3*N) = jerk_adim2;
jerk_adim_tot(3*N+1:4*N) = jerk_adim_dwell2;

% GEOMETRICAL
% Allocation of space
t_geom_tot = zeros(1,4*N);
pos_geom_tot = zeros(1,4*N);
vel_geom_tot = zeros(1,4*N);
acc_geom_tot = zeros(1,4*N);
jerk_geom_tot = zeros(1,4*N);

% alpha
t_geom_tot(1:N) = t_geom1;
t_geom_tot(N+1:2*N) = t_geom_dwell1+t_geom1(end);
t_geom_tot(2*N+1:3*N) = t_geom2+t_geom_tot(2*N);
t_geom_tot(3*N+1:4*N) = t_geom_dwell2+t_geom_tot(3*N);

% position
pos_geom_tot(1:N) = pos_geom1;
pos_geom_tot(N+1:2*N) = pos_geom1(end)+pos_geom_dwell1;
pos_geom_tot(2*N+1:3*N) = pos_geom2;
pos_geom_tot(3*N+1:4*N) = pos_geom_tot(3*N)+pos_geom_dwell2;

% velocity
vel_geom_tot(1:N) = vel_geom1;
vel_geom_tot(N+1:2*N) = vel_geom1(end)+vel_geom_dwell1;
vel_geom_tot(2*N+1:3*N) = vel_geom2;
vel_geom_tot(3*N+1:4*N) = vel_geom_tot(3*N)+vel_geom_dwell2;

% acceleration
acc_geom_tot(1:N) = acc_geom1;
acc_geom_tot(N+1:2*N) = acc_geom1(end)+acc_geom_dwell1;
acc_geom_tot(2*N+1:3*N) = acc_geom2;
acc_geom_tot(3*N+1:4*N) = acc_geom_tot(3*N)+acc_geom_dwell2;

% jerk
jerk_geom_tot(1:N) = jerk_geom1;
jerk_geom_tot(N+1:2*N) = jerk_geom_dwell1;
jerk_geom_tot(2*N+1:3*N) = jerk_geom2;
jerk_geom_tot(3*N+1:4*N) = jerk_geom_dwell2;

% conversion for negative cam (X is inverted)
pos_geom_tot=-pos_geom_tot;
vel_geom_tot=-vel_geom_tot;
acc_geom_tot=-acc_geom_tot;
jerk_geom_tot=-jerk_geom_tot;
pos_adim_tot=-pos_adim_tot;
vel_adim_tot=-vel_adim_tot;
acc_adim_tot=-acc_adim_tot;
jerk_adim_tot=-jerk_adim_tot;

% Important parameters for performance evaluation
ca_positive = max(acc_adim_tot);
ca_negative = min(acc_adim_tot);
cv = max(vel_adim_tot);
ck = max(abs(vel_adim_tot.*acc_adim_tot));

% Plots

figure
plot(t_geom_tot,pos_geom_tot)
grid on
hold on
plot(100,-160,'or')
hold on
plot(200,-160,'or')
xlim([0,t_geom_tot(end)])
title('Geometrical Position')
xlabel('\alpha [°]')
ylabel('y [mm]')

figure
plot(t_geom_tot,vel_geom_tot)
grid on
xlim([0,t_geom_tot(end)])
title('Geometrical Velocity')
xlabel('\alpha [°]')
ylabel(['y' char(39) '[mm/°]'])

figure
plot(t_geom_tot,acc_geom_tot)
grid on
xlim([0,t_geom_tot(end)])
title('Geometrical Acceleration')
xlabel('\alpha [°]')
ylabel(['y' char(39) char(39) '[mm/°^2]'])

figure
plot(t_geom_tot,jerk_geom_tot)
grid on
xlim([0,t_geom_tot(end)])
title('Geometrical Jerk')
xlabel('\alpha [°]')
ylabel(['y' char(39) char(39) char(39) '[mm/°^3]'])

figure
plot(t_adim_tot,pos_adim_tot)
grid on
hold on
plot(100/360,-160/180,'or')
hold on
plot(200/360,-160/180,'or')
xlim([0,t_adim_tot(end)])
title('Adimensional Position')
xlabel('\xi [-]')
ylabel('s(\xi) [-]')

figure
plot(t_adim_tot,vel_adim_tot)
grid on
xlim([0,t_adim_tot(end)])
title('Adimensional Velocity')
xlabel('\xi [-]')
ylabel('v(\xi) [-]')

figure
plot(t_adim_tot,acc_adim_tot)
grid on
xlim([0,t_adim_tot(end)])
title('Adimensional Acceleration')
xlabel('\xi [-]')
ylabel('a(\xi) [-]')

figure
plot(t_adim_tot,jerk_adim_tot)
grid on
xlim([0,t_adim_tot(end)])
title('Adimensional Jerk')
xlabel('\xi [-]')
ylabel('j(\xi) [-]')

figure
plot(t_adim_tot,vel_adim_tot.*acc_adim_tot)
grid on
title('ADIMENSIONAL V*A')
xlabel('\xi [-]')
ylabel('V*A [-]')
xlim([0,t_adim_tot(end)])

% conversion in radiants needed for cam synthesis
t_geom_tot = t_geom_tot*pi/180;
vel_geom_tot = vel_geom_tot*180/pi;
acc_geom_tot = acc_geom_tot*(180/pi)^2;
jerk_geom_tot = jerk_geom_tot*(180/pi)^3;

%% CAM SYNTHESIS %%%%%%%%%%%%%%%%%%%%%%%%%
% Here are y, y' and y'' needed for synthesis of cam
data.pos=pos_geom_tot;
data.vel=vel_geom_tot;
data.acc=acc_geom_tot;
alpha=t_geom_tot;

% Choose if you want eccentricity or not
eccentricityDesired = false;
e= 15;
if ~eccentricityDesired
    e = 0;
end

% First guess base radius for both cases with and without eccentricity
Rb0_e = h;
Rb0 = h;

% Roller radius choice
Rr = 20;

% Choose the limits for the pressure angle
theta_max_pos = 30;
theta_max_neg = 30;

% Loop for computation of base radius so that:
% 1. No Undercut
% 2. Max pressure angle limits are fulfilled

undercut_in=1;
undercut_out=1;
pressureAngleBig = true;
Undercut = true;
firstIteration = true;

while Undercut || pressureAngleBig
    if eccentricityDesired
        % Here important parameters for design
        if ~firstIteration
            Rb0_e = Rb0_e+1;
        end
        s0=Rb0_e; %INITIAL POSITION
        
        [xbaseradius,ybaseradius]=pol2cart(alpha,s0*ones(1,4*N)); %I just wanted to plot also the base radius
        
        
        % Polar coordinates of the pitch curve
        
        r0=sqrt((s0+data.pos).^2+e^2);
        phi0=alpha+atan(e./(s0+data.pos));
        
        [xpitch,ypitch]=pol2cart(phi0,r0);
        
        
        % pressure angle
        
        theta=atan((data.vel-e)./(s0+data.pos));
        
        % Polar coordinates of the cam profile
        
        r_in=sqrt((s0+data.pos-Rr*cos(theta)).^2+(e+Rr*sin(theta)).^2);
        phi_in=alpha+atan((e+Rr*sin(theta))./(s0+data.pos-Rr*cos(theta)));
        
        [xcam_in,ycam_in]=pol2cart(phi_in,r_in);
        
        r_out=sqrt((s0+data.pos+Rr*cos(theta)).^2+(e-Rr*sin(theta)).^2);
        phi_out=alpha+atan((e-Rr*sin(theta))./(s0+data.pos+Rr*cos(theta)));
        
        [xcam_out,ycam_out]=pol2cart(phi_out,r_out);
        
        
        % Curvature radius of the cam profile
        
        rho0=(((s0+data.pos).^2+(e-data.vel).^2).^(3/2))./((s0+data.pos).^2+...
            (e-data.vel).*(e-2*data.vel)-(s0+data.pos).*data.acc);
        rho_in=rho0-Rr;
        product=rho_in.*rho0;
        
        rho_out=rho0+Rr;
        product2=rho_out.*rho0;
        
    else
        % Here important parameters for design
        if ~firstIteration
            Rb0 = Rb0 + 1;
        end
        s0=Rb0; %INITIAL POSITION
        
        [xbaseradius,ybaseradius]=pol2cart(alpha,s0*ones(1,4*N)); %I just wanted to plot also the base radius
        
        
        % Polar coordinates of the pitch curve
        
        r0=s0+data.pos;
        phi0=alpha;
        
        [xpitch,ypitch]=pol2cart(phi0,r0);
        
        
        % pressure angle
        
        theta=atan(data.vel./r0);
        
        % Polar coordinates of the cam profile
        
        r_in=sqrt(Rr^2+r0.^2-2*Rr*r0.*cos(theta));
        phi_in=alpha+asin(Rr*sin(theta)./r_in);
        
        [xcam_in,ycam_in]=pol2cart(phi_in,r_in);
        
        r_out=sqrt(Rr^2+r0.^2+2*Rr*r0.*cos(theta));
        phi_out=alpha-theta+asin(r0.*sin(theta)./r_out);
        
        [xcam_out,ycam_out]=pol2cart(phi_out,r_out);
        
        
        % Curvature radius of the cam profile
        
        rho0=(data.vel.^2+r0.^2).^(3/2)./(r0.^2-r0.*data.acc+2*data.vel.^2);
        rho_in=rho0-Rr;
        product=rho_in.*rho0;
        
        rho_out=rho0+Rr;
        product2=rho_out.*rho0;
    end
    
    % Undercut assessment
    
    undercut_in=0;
    for i=1:length(product)
        if product(i)<0
            undercut_in=1;
        end
    end
    
    undercut_out=0;
    for i=1:length(product2)
        if product2(i)<0
            undercut_out=1;
        end
    end
    
    if undercut_in==0 && undercut_out==0
        Undercut = false;
    end
    
    if max(theta)<theta_max_pos*pi/180 && min(theta)>-theta_max_neg*pi/180
        pressureAngleBig = false;
    end
    firstIteration = false;
end

% Plot of the pitch and cam profile

figure
plot(xpitch,ypitch,'Displayname','Pitch')
hold on
plot(xcam_in,ycam_in,'Displayname','Cam in')
hold on
plot(xcam_out,ycam_out,'Displayname','Cam out')
hold on
plot(xbaseradius,ybaseradius,'Displayname','Base Radius')
grid on
legend
xlabel('X [mm]')
ylabel('Y [mm]')
title('Profiles')
xlim ([-1.25*max(r_out) 1.25*max(r_out)]) %these 3 last rows of code are to center the plot and keep the same axis proportion
ylim ([-1.25*max(r_out) 1.25*max(r_out)])
daspect([1 1 1])


% Plot of the pressure angle

figure
plot(t_geom_tot*180/pi,theta*180/pi)
grid on
xlabel('\alpha [deg]')
ylabel('theta [deg]')
title('Pressure angle')

% Plots of rho_in, rho_out and rho_0
figure
plot(t_geom_tot*180/pi,rho0)
grid on
title('Curvature radius pitch profile')
xlabel('\alpha [deg]')
ylabel('\rho_0 [mm]')

figure
plot(t_geom_tot*180/pi,rho_in)
grid on
title('Curvature radius internal profile')
xlabel('\alpha [deg]')
ylabel('\rho_{in} [mm]')

figure
plot(t_geom_tot*180/pi,rho_out)
grid on
title('Curvature radius external profile')
xlabel('\alpha [deg]')
ylabel('\rho_{out} [mm]')


% Plots of the product rho_0*rho for internal and external profile

figure
plot(t_geom_tot*180/pi,product)
grid on
if undercut_in==0
    title('\rho e \rho_0 concordi: NO UNDERCUT IN')
else
    title('\rho e \rho_0 discordi: YES UNDERCUT IN')
end
legend('\rho\rho_0')
xlabel('\alpha [deg]')
ylabel('\rho\rho_0 [mm]')


figure
plot(t_geom_tot*180/pi,product2)
grid on
if undercut_out==0
    title('\rho e \rho_0 concordi: NO UNDERCUT OUT')
else
    title('\rho e \rho_0 discordi: YES UNDERCUT OUT')
end
legend('\rho\rho_0')
xlabel('\alpha [deg]')
ylabel('\rho\rho_0 [mm]')



%% Profiles to be exported
internalCam = zeros(4*N,3);
externalCam = zeros(4*N,3);
pitchProfile = zeros(4*N,3);

internalCam(:,1) = xcam_in;
internalCam(:,2) = ycam_in;
internalCamTab = table(internalCam);

externalCam(:,1) = xcam_out;
externalCam(:,2) = ycam_out;
externalCamTab = table(externalCam);

pitchProfile(:,1) = xpitch;
pitchProfile(:,2) = ypitch;
pitchProfileTab = table(pitchProfile);

writetable(internalCamTab,'internalCam.txt');
writetable(externalCamTab,'externalCam.txt');
writetable(pitchProfileTab,'pitchProfileFirstCam.txt');



%% Saving motionlaw for design second cam
% we save the motionlaw only of the rise for the design of the second cam

% we invert the sign because we want the other code considered the correct
% sign

data_x.pos=-pos_geom_tot(1:N);
data_x.vel=-vel_geom_tot(1:N);
data_x.acc=-acc_geom_tot(1:N);
data_x.jerk=-jerk_geom_tot(1:N);

save('data_x.mat','data_x');

%% Saving curvature radius for contact Pressure computation
curvature_radii_firstCam.rho_in = rho_in;
curvature_radii_firstCam.rho_out = rho_out;
curvature_radii_firstCam.rho0 = rho0;

save('curvature_radii_firstCam.mat','curvature_radii_firstCam');
