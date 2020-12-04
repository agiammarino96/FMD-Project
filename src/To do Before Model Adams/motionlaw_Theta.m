close all
clc
clear all
MotionLaw = {'Dwell' 'CIC' 'Dwell' 'CIC' 'Dwell'};

LiftH = [0 +90 0 -90 0];
Abscissa = [100 30 40 30 160]; % duration of each phase

N = 400; % total number of discretization intervals
abscissa = linspace (0 ,sum( Abscissa ) ,N); % vector with N elements

delta=[0.125, 0, 0.375, 0, 0.375, 0, 0.125]; 
data = motionlaw(MotionLaw , LiftH , Abscissa , N, delta);		
figure
subplot(4,1,1);
plot(abscissa,data.pos,'Displayname','CIC');
ylabel('Theta.[ deg ]');grid ;
xlim([abscissa(1),abscissa(end)]);
title ('MotionLaw');
legend
subplot(4,1,2); plot(abscissa,data.vel); 
ylabel('Vel.[ deg/ab ]'); grid ;
xlim([abscissa(1),abscissa(end)]);
subplot(4,1,3); plot(abscissa,data.acc);
ylabel('Acc.[ deg/ab\^2 ]');  grid ;	
xlabel('Abscissa');
xlim([abscissa(1),abscissa(end)]);
subplot(4,1,4); plot(abscissa,data.jerk);
ylabel('Jerk.[ deg/ab\^3 ]');  grid ;	
xlabel('Abscissa');
xlim([abscissa(1),abscissa(end)]);

%% Saving motionlaw
% we want to save only the rise for the second cam design
MotionLaw_rise = {'Dwell' 'CIC'};

LiftH_rise = [0 +90];
Abscissa_rise = [100 30]; % duration of each phase


N_rise = N/4; % it is important to have the same number in the motionlaw of x
abscissa = linspace (0 ,sum( Abscissa_rise ) ,N_rise); % vector with N elements

delta=[0.125, 0, 0.375, 0, 0.375, 0, 0.125]; 
data_theta = motionlaw(MotionLaw_rise , LiftH_rise , Abscissa_rise , N_rise, delta);
data_theta.pos=data_theta.pos*pi/180;
data_theta.vel=data_theta.vel;
data_theta.acc=data_theta.acc*(180/pi);
data_theta.jerk=data_theta.jerk*(180/pi)^2;


save('data_theta.mat','data_theta')
