%% Computation contact Pressure (SCRIPT NOT COMPLETED)
clear all
close all
clc

% PAY ATTENTION: the contact forces you import from Adams will have a
% certain number of samples and will be referred to a certain time interval
% (corresponding to one or more rotations of the cam). The code should be
% written in order to adapt the curvature radii vectors and the contact
% force vectors since they are used in the same formula in order to compute
% the contact pressure. Suggestion: adapt the curvature radii to the Adams
% files since you can generate them easily from the script changing the
% number of points used.

%% These files are generated from the scripts for cam design
load('curvature_radii_secondCam.mat');
load('curvature_radii_firstCam.mat');

rho_in_firstCam = curvature_radii_firstCam.rho_in;
rho_out_firstCam = curvature_radii_firstCam.rho_out;
rho_out_secondCam = curvature_radii_secondCam.rho_out;
rho_in_secondCam = curvature_radii_secondCam.rho_in;

N = length(rho_out_firstCam);
time_cam = linspace(0,2,N); % 2 seconds needed for a turn

%% Import file from Adams containing contact forces
A = importdata('ExternalContactForceCam1.txt');
B = importdata('ExternalContactForceCam2.txt');
C = importdata('InternalContactForceCam1.txt');
% D = importdata('InternalContactForceCam2.txt');
E = importdata('FrictionForceExternalContact.txt');
G = importdata('FrictionForceInternalContact.txt');
time = A.data(:,1);
externalContactForceCam1 = A.data(:,2);
externalContactForceCam2 = B.data(:,2);
internalContactForceCam1 = C.data(:,1);
% internalContactForceCam2 = D.data(:,2);
totalFrictionExternalContact = E.data(:,1)+E.data(:,2)+E.data(:,3)+E.data(:,4);
totalFrictionInternalContact = G.data(:,1)+G.data(:,2)+G.data(:,3)+G.data(:,4)+G.data(:,5);

figure
plot(time,totalFrictionExternalContact)
hold on
grid on
plot(time,externalContactForceCam1*0.3)
title('Slipping check External Contact Cam1')
legend('Friction Force','Max Friction allowable')

figure
plot(time,totalFrictionInternalContact)
hold on
grid on
plot(time,internalContactForceCam1*0.3)
title('Slipping check Internal Contact Cam1')
legend('Friction Force','Max Friction allowable')


%% Contact pressure computation First Cam
E = 206000; % Young Modulus [MPa]
b1 = 20; % Roller thickness [mm]
Rr1 = 20; % Roller radius [mm]
Rr1_vector = ones(1,N)*Rr1;

i = find(time==2.0); % Initial time Adams simulation
externalContactForceCam1OneTurn = externalContactForceCam1(i:i+N-1)'; % contact force corresponding to one turn

P_contact1 = sqrt((0.175*E*externalContactForceCam1OneTurn)./(b1*Rr1_vector).*(1+Rr1_vector./rho_out_firstCam)); % [MPa]
Pmax_cam1 = max(P_contact1);
disp(['Pmax first cam = ', num2str(Pmax_cam1), ' MPa'])

figure
plot(time_cam,P_contact1)
grid on
title('Contact Pressure External Profile first cam')
xlabel('time [s]')
ylabel('P_{contact,1} [MPa]')

%% Contact pressure computation Second Cam
E = 206000; % Young Modulus [MPa]
b2 = 20; % Roller thickness [mm]
Rr2 = 10; % Roller radius [mm]
Rr2_vector = ones(1,N)*Rr2;

externalContactForceCam2OneTurn = externalContactForceCam2(i:i+N-1)'; % contact force corresponding to one turn

P_contact2 = sqrt((0.175*E*externalContactForceCam2OneTurn)./(b2*Rr2_vector).*(1+Rr2_vector./rho_out_firstCam)); % [MPa]
Pmax_cam2 = max(P_contact2);
disp(['Pmax second cam = ', num2str(Pmax_cam2), ' MPa'])

figure
plot(time_cam,P_contact2)
grid on
title('Contact Pressure External Profile second cam')
xlabel('time [s]')
ylabel('P_{contact,2} [MPa]')

%% MotorTorque Data
% We should take only the torque of the steady state part, so excluding the
% initial transient
T = importdata('MotorTorque.txt');
time = T.data(:,1);
Torque = T.data(:,2);
i = find(time==0.2); % In this case transient was 0.1s
time = time(i:end);
Torque = Torque(i:end)*10^(-3);
Torque_rms = rms(Torque);
Torque_max = max(abs(Torque));
omega = 3.14; % [rad/s] Motor speed

figure
plot(time,Torque)
grid on
hold on
plot([time(1) time(end)],[Torque_rms Torque_rms])
title('Motor Torque')
xlabel('time [s]')
ylabel('Torque [N*m]')
legend('Torque Steady State','Torque RMS')
xlim([time(1) time(end)])

% Motor parameters
Ratio_Cmax_Cn_motor = 2.6;
C_n_motor = 1.26; % [N*m] Nominal torque motor
C_max_motor = Ratio_Cmax_Cn_motor*C_n_motor; % [N*m] Maximum torque motor
omega_n_motor = 2815/60; % [rad/s] Nominal speed motor

% Transmission choice and transformation motor parameters
tau = 1/8;
C_max_load = C_max_motor/tau;
C_n_load = C_n_motor/tau;
omega_n_load = omega_n_motor*tau;

% Check
if Torque_max<C_max_load
    disp('Max Torque Check ok')
    disp([num2str(Torque_max), ' < ', num2str(C_max_load)])
else
    disp('MAX TORQUE CHECK NOT SUCCEEDED')
    disp([num2str(Torque_max), ' > ', num2str(C_max_load)])
end

if Torque_rms<C_n_load
    disp('RMS Torque Check ok')
    disp([num2str(Torque_rms), ' < ', num2str(C_n_load)])
else
    disp('RMS TORQUE CHECK NOT SUCCEEDED')
    disp([num2str(Torque_rms), ' > ', num2str(C_n_load)])
end

if omega<omega_n_load
    disp('Speed Check ok')
    disp([num2str(omega), ' < ', num2str(omega_n_load)])
else
    disp('SPEED CHECK NOT SUCCEEDED')
    disp([num2str(omega), ' > ', num2str(omega_n_load)])
end
