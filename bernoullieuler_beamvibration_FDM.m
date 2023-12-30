% Bernoulli-Euler Beam Vibration with FDM

% Parameters

L=41.4*10^-2; % m  % The length of the beam
wd=14.88*10^-3; % m  % The width of the beam
t=1*10^-3; % m    % The thickness of the beam 
E=70*10^9; % N/m^2   % The modulus of elasticity of the beam
ro=5186.5; % kg/m^2 The linear mass density of the beam
I=((wd*(t^3))/12); % kgm^2 The moment of inertia of the beam
A=wd*t; % m^2 The cross-sectional area of the beam
lamda=A*ro ; 
ksi=0.2; % The damping ratio of the beam
P=5*9.81; % Tip payload
nt=10^5; % length of time domain
nx=15; % length of space domain
tf=15; % time of simulation

% Compute the mesh spacing and time step

dx=L/nx; % spacing step 
dt=tf/nt; % time step

% Drawing coordinate

x=linspace(0,L,nx+1); % create spaceline for drawing
t=linspace(0,tf,nt+1); % create timeline for drawing

% Create memory to save data contains

w=zeros(nx+1,nt+1); % vibrations at position x for time t
d=zeros(nt+1,1); % boundary disturbance
f=zeros(nx+1,nt+1); % disturbed disturbance

k1=1.875104069/L;
k2=4.694091133/L;
k3=7.854757438/L;
k4=10.99554073/L; 
k5=14.13716839/L;
k6=17.27875953/L;

w1=sqrt((((k1)^4)*E*I)/(lamda));
w2=sqrt((((k2)^4)*E*I)/(lamda));
w3=sqrt((((k3)^4)*E*I)/(lamda));
w4=sqrt((((k4)^4)*E*I)/(lamda));
w5=sqrt((((k5)^4)*E*I)/(lamda));
w6=sqrt((((k6)^4)*E*I)/(lamda));

ra=(2*ksi*w1*ro*A);

% Boundary disturbance

for j=1:nt+1
d(j)=0.1+0.1*sin(pi*(j-1)*dt)+0.1*sin(2*pi*(j-1)*dt)+0.1*sin(3*pi*(j-1)*dt);
end

% Disturbed disturbance

for i=1:nx+1
for j=1:nt+1
f(i,j)=(1+sin(0.1*pi*(i-1)*dx*(j-1)*dt)+sin(0.2*pi*(i-1)*dx*(j-1)*dt)+sin(0.3*pi*(i-1)*dx*(j-1)*dt))*(i-1)*dx/(20*L);
end
end

% Initial Condition's

for i=1:nx+1
w(i,1)=-(P*L/(E*I))*(((x(i).^2)./2)-((x(i).^3)./(6*L)));
w(i,2)=w(i,1);
end

alpha=(dt*ra/(ro*A));

T=6; % Tension
Ms=5; % Mass of the tip payload

w(1,:)=0;
w(2,:)=w(1,:);

% Main loop

for j=3:nt+1     
for i=3:nx-1
wxxxx=(w(i+2,j-1)-4*w(i+1,j-1)+6*w(i,j-1)-4*w(i-1,j-1)+w(i-2,j-1))/dx^4;
w(i,j)=((2+alpha)*w(i,j-1)-w(i,j-2)+(-E*I*wxxxx+f(i,j))*(dt^2/(ro*A)))/(alpha+1);
end

% Boundary Condition's
wxL=(w(nx+1,j-1)-w(nx,j-1))/dx;
wxxxL=(w(nx+1,j-1)-3*w(nx,j-1)+3*w(nx-1,j-1)-w(nx-2,j-1))/dx^3;
w(nx+1,j)=2*w(nx+1,j-1)-w(nx+1,j-2)+(d(j-1)+E*I*wxxxL-T*wxL)*dt^2/Ms;
w(nx,j)=(w(nx+1,j)+w(nx-1,j))/2; 
end

% Make a draw

figure (1)
mesh(x,t,w');
title('Displacement of beam without control');
ylabel('Time(s)','Fontsize',12);
xlabel('x(m)','Fontsize',12);
zlabel('w(x,t)(m)','Fontsize',12);
view([60 45]);

% Animation
figure;
for j = 1:length(t)
    plot(x, w(:,j));
    xlabel('x (m)');
    ylabel('Displacement (m)');
    title(sprintf('Beam Displacement at Time = %.2f s', t(j)));
    axis([0 L min(w(:)) max(w(:))]);
    drawnow;
    pause(0.01); % pause for a short time to slow down the animation
end