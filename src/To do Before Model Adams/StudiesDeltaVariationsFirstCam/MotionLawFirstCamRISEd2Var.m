close all
clear all
clc
%% Useful vectors and parameters First part
% Choose d1, d2, d3 and d4; compute d5, d6 and d7
M = 4;
N = 1000; % Number of samples
h = 180; % Total rise motion law
Ab = 130; % Cam angle correspondent to rise
d1_var = zeros(1,M);
ca_negative1_var = zeros(1,M);
ca_positive1_var = zeros(1,M);
ck1_var = zeros(1,M);
cv1_var = zeros(1,M);
cj1_var = zeros(1,M);
pos_adim1 = zeros (M,N);
vel_adim1 = zeros (M,N);
acc_adim1 = zeros (M,N);
jerk_adim1 = zeros(M,N);
pos_geom1 = zeros (M,N);
vel_geom1 = zeros (M,N);
acc_geom1 = zeros (M,N);
jerk_geom1 = zeros(M,N);
t_geom1 = zeros(M,N);
t_adim1 = zeros(M,N);

for j=1:M
    d2_var(j)=-0.02+0.02*j;
    d = zeros(1,7);
    d(1) = 0.22;
    d(2) = d2_var(j);
    d(3) = 0.22;
    d(4) = 0.0;
    d(5) = 1-(30/130+sum(d(1:4))); % chosen so that the sum(d(1:5)) is so that alpha_5=100° (d6+d7=30/130)
    
    t_geom1(j,:) = linspace (0, Ab, N); % Cam angle discretization
    t_adim1(j,:) = linspace(0,1,N); % Discretization for adimensional motion law
    %% Computation of A, B, d2 and d6 so that
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
    
    %% Values of A, B and d computed converted to double
    A_solution = double(A);
    B_solution = double(B);
    d6_solution = double(d6);
    d7_solution = double(d7);
    
    %% Choice of the solution to use
    A = A_solution(2);
    B = B_solution(2);
    d(6) = d6_solution(2);
    d(7) = d7_solution(2);
    
    %% Computation of motion law
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
    
    for i=1:length(d)
        a(i)=sum(d(1:i));
    end
    
    for k=1:N
        epsilon = t_geom1(j,k)/Ab;
        
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
        pos_geom1(j,k) = out.s*h;
        vel_geom1(j,k) = out.v*h/Ab;
        acc_geom1(j,k) = out.a*h/Ab^2;
        jerk_geom1(j,k) = out.j*h/Ab^3;
        pos_adim1(j,k) = out.s;
        vel_adim1(j,k) = out.v;
        acc_adim1(j,k) = out.a;
        jerk_adim1(j,k) = out.j;
    end
    
    
    %% Computation ca- ca+ ck and cv
    
    ca_positive1 = max(acc_adim1(j,:));
    ca_negative1 = min(acc_adim1(j,:));
    [cv1, index] = max(vel_adim1(j,:));
    ck1 = max(abs(vel_adim1(j,:).*acc_adim1(j,:)));
    cj1 = abs((jerk_adim1(j,index-1)+jerk_adim1(j,index))/2);
    
    ca_positive1_var(j) = ca_positive1;
    ca_negative1_var(j) = ca_negative1;
    cv1_var(j) = cv1;
    ck1_var(j) = ck1;
    cj1_var(j) = cj1;
    
    % %% Plots Geometrical laws
    % figure(1)
    % subplot(4,1,1)
    % plot(t_geom1,jerk_geom1)
    % grid on
    % title('GEOMETRICAL LAWS RISE')
    % xlabel('\alpha [°]')
    % ylabel('Jerk [mm/°^3]')
    % xlim([0,Ab])
    %
    % subplot(4,1,2)
    % plot(t_geom1,acc_geom1)
    % grid on
    % xlabel('\alpha [°]')
    % ylabel('Acc [mm/°^2]')
    % xlim([0,Ab])
    %
    % subplot(4,1,3)
    % plot(t_geom1,vel_geom1)
    % grid on
    % xlabel('\alpha [°]')
    % ylabel('Vel [mm/°]')
    % xlim([0,Ab])
    %
    % subplot(4,1,4)
    % plot(t_geom1,pos_geom1)
    % hold on
    % plot(100,160,'or')
    % grid on
    % xlabel('\alpha [°]')
    % ylabel('Pos [mm]')
    % xlim([0,Ab])
    
    %% Plots Adimensional laws
    
%     figure(2)
%     subplot(4,1,1)
%     plot(t_adim1,jerk_adim1)
%     grid on
%     hold on
%     title('ADIMENSIONAL LAWS RISE')
%     xlabel('\xi [-]')
%     ylabel('Jerk [-]')
%     xlim([0,1])
%     
%     subplot(4,1,2)
%     plot(t_adim1,acc_adim1)
%     grid on
%     hold on
%     xlabel('\xi [-]')
%     ylabel('Acc [-]')
%     xlim([0,1])
%     
%     subplot(4,1,3)
%     plot(t_adim1,vel_adim1)
%     grid on
%     hold on
%     xlabel('\xi [-]')
%     ylabel('Vel [-]')
%     xlim([0,1])
%     
%     subplot(4,1,4)
%     plot(t_adim1,pos_adim1)
%     grid on
%     hold on
%     xlabel('\xi [-]')
%     ylabel('Pos [-]')
%     xlim([0,1])
%     
    
end

%% Plots
    figure
    for j=1:M
        plot(t_adim1(j,:),vel_adim1(j,:).*acc_adim1(j,:))
        grid on
        hold on
    end
    title('ADIMENSIONAL V*A RISE')
    xlabel('\xi [-]')
    ylabel('V*A [-]')
    xlim([0,1])
    legend('1','2','3','4')
    
    figure
    subplot(2,1,1)
    for j=1:M
        plot(t_adim1(j,:),jerk_adim1(j,:))
        grid on
        hold on
    end
    title('ADIMENSIONAL ACC AND JERK')
    xlabel('\xi [-]')
    ylabel('Jerk [-]')
    xlim([0,1])
    legend('1','2','3','4')
    
    subplot(2,1,2)
    for j=1:M
        plot(t_adim1(j,:),acc_adim1(j,:))
        grid on
        hold on
    end
    xlabel('\xi [-]')
    ylabel('Acc [-]')
    xlim([0,1])
    legend('1','2','3','4')
    
    figure
    for j=1:M
        plot(t_adim1(j,:),vel_adim1(j,:))
        grid on
        hold on
    end
    title('ADIMENSIONAL VEL')
    xlabel('\xi [-]')
    ylabel('Vel [-]')
    xlim([0,1])
    legend('1','2','3','4')

    figure
    plot(d2_var,ca_negative1_var,'o-')
    hold on
    plot(d2_var,ca_positive1_var,'*-')
    title('Trends ca+ and ca-')
    xlabel('d2 [mm]')

    figure
    plot(d2_var,ck1_var,'o-')
    title('Trend ck')
    xlabel('d2 [mm]')
    
    figure
    plot(d2_var,cv1_var,'o-')
    title('Trend cv')
    xlabel('d2 [mm]')
    
    figure
    plot(d2_var,cj1_var,'o-')
    title('Trend cj')
    xlabel('d2 [mm]')
    
    