function [out] = motionlaw (MotionLaw , h , Ab, N, delta)

if length(MotionLaw) ~=length(h) || length(h) ~= length(Ab)
    error ('MotionLaw , h and Ab arrays must have the same size' );
end

t = linspace (0,sum(Ab),N);

% Dimensional case (angle domain)
out . pos = zeros (1 ,N);
out . vel = zeros (1 ,N);
out . acc = zeros (1 ,N);
out . jerk = zeros (1 ,N);
out . specpow = zeros (1 ,N);
% Adimensional case
out . pos_adim = zeros (1 ,N);
out . vel_adim = zeros (1 ,N);
out . acc_adim = zeros (1 ,N);
out. jerk_adim = zeros (1, N);
out . specpow_adim = zeros (1 ,N);

PreviousCurInt = 1; % variable that stores the "CurInt" of previous iteration
shiftPosition = 0; % variable that corresponds to the position shift needed for continuity
shiftPositionAdim = 0;
for k = 1:N
    
    CurInt = find(cumsum(Ab)>=t(k),1,'first');
    
    % Check if CurInt has changed
    if PreviousCurInt == CurInt
        changedCurInt = false;
    else
        changedCurInt = true;
    end
    
    % Update Previous CurInt for next iteration
    PreviousCurInt = CurInt;
    
    
    epsilon = (t(k) - sum(Ab(1:(CurInt-1)))) / Ab(CurInt);
    
   
    switch char (MotionLaw(CurInt))
        case {'ACS'}
            AdimLaw = ACS(epsilon);
        case {'TVP'}
            AdimLaw = TVP(epsilon);
        case {'CUB'}
            AdimLaw = Cubic(epsilon);
        case {'CIC'}
            AdimLaw = Cycloidal(epsilon);
        case {'Dwell'}
            AdimLaw = Dwell(epsilon);
        case {'ModTVP'}
            AdimLaw = ModTVP(epsilon,delta);
        case {'ModTVPi'} % this  is done to have a symmetric motionlaw with an asymmetric delta
            Delta=fliplr(delta);
            AdimLaw = ModTVP(epsilon,Delta);
        otherwise
            error('s  motion law not recognized',char(MotionLaw(CurInt)));
    end
    
    % If CurInt has changed, define the shiftPosition as the last position
    % of the previous elementary law
    if changedCurInt
        shiftPosition = out.pos(k-1);
        shiftPositionAdim = out.pos_adim(k-1);
    end
    
   % Position is shifted for continuity
    out.pos(k) = AdimLaw.s*h(CurInt) + shiftPosition;
    out.vel (k) =AdimLaw.v*h(CurInt)/Ab(CurInt);
    out.acc(k) =AdimLaw.a*h(CurInt)/Ab(CurInt)^2;
    out.jerk(k) = AdimLaw.j*h(CurInt)/Ab(CurInt)^3;
    out.specpow(k) = out.vel (k)* out.acc(k);
    
    out.pos_adim(k) = AdimLaw.s + shiftPositionAdim;
    out.vel_adim(k) =AdimLaw.v;
    out.acc_adim(k) =AdimLaw.a;
    out.jerk_adim(k) = AdimLaw.j;
    out.specpow_adim(k) = out.vel_adim(k)* out.acc_adim(k);
end
end



function [out] = TVP(epsilon)

ev = 1/3;
cv = 1/(1-ev);
ca = 1/ev/(1-ev);


if 0 <= epsilon && epsilon <=ev
    out.j = 0;
    out.a =  ca;
    out.v =  ca*epsilon;
    out.s =  1/2*ca*epsilon^2;
    out.j = 0;
elseif ev < epsilon && epsilon <=(1-ev)
    out.j = 0;
    out.a = 0;
    out.v = ca*ev;
    out.s = ca*ev*epsilon-1/2*ca*ev^2;
else
    out.j = 0;
    out.a = -ca;
    out. v = ca*(1-epsilon);
    out. s = ca*(epsilon-epsilon^2/2+ev-ev^2-1/2);
end
end

function [out] = Cubic(epsilon)
delta = 2/3;
cv = 1/delta;
ca = 4*cv;
out.j = -2*ca;
out.a = ca*(1-2*epsilon);
out. v = ca*epsilon*(1-epsilon);
out. s = epsilon*(3*epsilon-2*epsilon^2);
end

function [out] = Dwell(epsilon)
out.j = 0;
out.a = 0;
out. v = 0;
out. s = 0;
end

function [out] = Cycloidal(epsilon)
    delta = 0.5;
    cv = 1/delta;
    ca = 2*pi;
    x = 2*pi*epsilon;
    
    out.j = ca*2*pi*cos(2*pi*epsilon);
    out.a =  ca*sin(x);
    out.v =  1 - cos(x);
    out.s =  epsilon - (1/ca)*sin(x);
end

function [out] = ACS(epsilon)

 ev = 1/2;
 cv = 2;
 % ca = max(2/ev,2/(1-ev));
 
if  epsilon <= ev	
    out.j = 0;
    out.a =  cv/ev;
    out.v =  (cv/ev)*epsilon;
    out.s =  1/2*(cv/ev)*epsilon^2;
else
     out.j = 0;
     out.a = -cv/(1-ev);
     out. v = cv/(1-ev)*(1-epsilon);
     out. s = cv/(1-ev)*(epsilon-epsilon^2/2-ev/2);
end
end

function [out] = ModTVP(epsilon,delta)

% remember to enter in the input a vector delta composed of 7 intervals

a11=2*delta(1)/pi+delta(2)+2*delta(3)/pi;
a12=-2*delta(5)/pi-delta(6)-2*delta(7)/pi;
a21=2*delta(1)/pi*((pi-2)/pi*delta(1)+delta(2)/2)+(2*delta(1)/pi+delta(2))*(delta(2)/2+(pi-2)/pi*delta(3))+(2*delta(1)/pi+delta(2)+2*delta(3)/pi)*(2*delta(3)/pi+delta(4)+2*delta(5)/pi);
a22=(2*delta(7)/pi+delta(6))*((pi-2)/pi*delta(5)+delta(6)/2)+2*delta(7)/pi*(delta(6)/2+(pi-2)/pi*delta(7));
A=-a12/(a11*a22-a12*a21);
B=a11/(a11*a22-a12*a21);
Epsilon=cumsum(delta);
v1 = A*2*delta(1)/pi;
s1 = A*2*delta(1)^2/pi*(1-2/pi);
v2 = v1+A*delta(2);
s2 = s1+A*delta(2)*(2*delta(1)/pi+delta(2)/2);
v3 = v2+A*2*delta(3)/pi;
s3 = s2+A*delta(3)*(2*delta(1)/pi+delta(2)+4*delta(3)/(pi^2));
v4 = v3;
s4 = s3+A*delta(4)*(2*delta(1)/pi+delta(2)+2*delta(3)/pi);
v5 = v4-B*2*delta(5)/pi;
s5 = s4+A*delta(5)*(2*delta(1)/pi+delta(2)+2*delta(3)/pi)-B*2*delta(5)^2/pi*(1-2/pi);
v6 = v5-B*delta(6);
s6 = s5+v5*delta(6)-B*delta(6)^2/2;
v7 = v6-B*2*delta(7)/pi;
s7 = s6+v6*delta(7)-B*(2*delta(7)/pi)^2;

 

if 0 <= epsilon && epsilon <=Epsilon(1)	
    out.j = A*pi/(2*delta(1)) * cos(epsilon*pi/(2*delta(1)));
    out.a = A*sin(epsilon*pi/(2*delta(1)));
    out.v = A*2*delta(1)/pi*(1-cos(epsilon*pi/(2*delta(1))));
    out.s = A*2*delta(1)/pi*(epsilon-2*delta(1)/pi*sin(epsilon*pi/(2*delta(1))));
    
elseif Epsilon(1) < epsilon && epsilon <=Epsilon(2)
     out.j = 0;
     out.a = A;
     out.v = v1+A*(epsilon-Epsilon(1));
     out.s = s1+v1*(epsilon-Epsilon(1))+A*(epsilon-Epsilon(1))^2/2;
     
elseif Epsilon(2) < epsilon && epsilon <=Epsilon(3)
     out.j = -A*pi/(2*delta(3)) * sin((epsilon-Epsilon(2))*pi/(2*delta(3)));
     out.a = A*cos((epsilon-Epsilon(2))*pi/(2*delta(3)));
     out.v = v2+A*2*delta(3)/pi*sin((epsilon-Epsilon(2))*pi/(2*delta(3)));
     out.s = s2+v2*(epsilon-Epsilon(2))+A*(2*delta(3)/pi)^2*(1-cos((epsilon-Epsilon(2))*pi/(2*delta(3))));
     
elseif Epsilon(3) < epsilon && epsilon <=Epsilon(4)
     out.j = 0;
     out.a = 0;
     out.v = v3;
     out.s = s3+v3*(epsilon-Epsilon(3));
     
elseif Epsilon(4) < epsilon && epsilon <=Epsilon(5)
     out.j = -B*pi/(2*delta(5)) * cos((epsilon-Epsilon(4))*pi/(2*delta(5)));
     out.a = -B*sin((epsilon-Epsilon(4))*pi/(2*delta(5)));
     out.v = v4-B*2*delta(5)/pi*(1-cos((epsilon-Epsilon(4))*pi/(2*delta(5))));
     out.s = s4+v4*(epsilon-Epsilon(4))-B*2*delta(5)/pi*(epsilon-Epsilon(4)-2*delta(5)/pi*sin((epsilon-Epsilon(4))*pi/(2*delta(5))));
     
elseif Epsilon(5) < epsilon && epsilon <=Epsilon(6)
     out.j = 0;
     out.a = -B;
     out.v = v5-B*(epsilon-Epsilon(5));
     out.s = s5+v5*(epsilon-Epsilon(5))-B*(epsilon-Epsilon(5))^2/2;
     
else
     out.j = B*pi/(2*delta(7)) * sin((epsilon-Epsilon(6))*pi/(2*delta(7)));
     out.a = -B*cos((epsilon-Epsilon(6))*pi/(2*delta(7)));
     out.v = v6-B*2*delta(7)/pi*sin((epsilon-Epsilon(6))*pi/(2*delta(7)));
     out.s = s6+v6*(epsilon-Epsilon(6))-B*(2*delta(7)/pi)^2*(1-cos((epsilon-Epsilon(6))*pi/(2*delta(7))));
     
end
end

