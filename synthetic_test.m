% Numerical simulation of the reflection and transmission of an elastic wave 
% normal to a fluid-filled layer, and fluid estimation using the feq-domain
% technique from Groenenboom paper.

close all
clearvars -EXCEPT P_*
home

%% Problem definition
% meshing and constants
c0 = 2800;
dx = 0.4*c0*1E-8;
dt = 0.4E-8;
Fn = 1/(2*dt);

f0 = 1E6;
X = 0:dx:0.3;
T = (0:dt:1E-4)';

% homogeneous medium
C = ones(size(X))*c0*dt/dx;
c1 = 1250;
h = 2E-5; %1E-4;
x0 = 0.125;
% heterogeneous medium
C(X>=x0 & X<=x0+h) = c1/c0;

% initial conditions
P0 = zeros(size(X));
Q0 = zeros(size(X));
P1 = zeros(size(X));
Q1 = zeros(size(X));
P2 = zeros(size(X));

% boundary conditions
pinit = zeros(size(T));
w = ricker(f0,dt);
pinit(1:length(w)) = w;
P0(1) = pinit(1);

%% computation
% time n=2 with Lax-wendroff
for nn = 2
    P1(1)=pinit(nn);
    for j = 2:length(X)-1
        P1(j) = P0(j) - C(j)/2 * (Q0(j+1)-Q0(j-1)) + ...
            C(j)^2 /2 *(P0(j-1)-2*P0(j)+P0(j+1)) ;
        Q1(j) = Q0(j) - C(j)/2 * (P0(j+1)-P0(j-1)) + ...
            C(j)^2 /2 *(Q0(j-1)-2*Q0(j)+Q0(j+1)) ;
    end
end
figure

ii_end = find(X>=2*x0+h,1);
P_end = zeros(size(T));
P_end(1) = P0(ii_end);
P_end(2) = P1(ii_end);

% time n>2 with Leap-Frog
for nn = 3:length(T)
    P2(1) = pinit(nn);
    for j = 2:length(X)-1
        P2(j) = 2*P1(j)-P0(j)+C(j)^2 * (P1(j+1)+P1(j-1)-2*P1(j));
    end
    % plotting at every kk step
    kk = 20;
    if mod(nn,kk) == 0
        plot(X*1E2,P2,'r',[x0 x0]*1E2,[-1 1]*2,':k')
        title(['iteration number ' num2str(nn)])
        %axis([x(1) x(length(x)) -1.1 1.1])
        axis([[X(1) x0*2+h]*1E2 -1.1 1.1])
        xlabel('Position (cm)')
        ylabel('Pressure (a.u.')
        pause(0.01)
    end
    
    P_end(nn) = P2(ii_end);
    P0 = P1;
    P1 = P2;
    P2 = zeros(size(X));
end

%% check simulation at end of block
figure
plot(T*1E6,P_end)
axis([[T(1) T(end)]*1E6 [-1 1]])
xlabel('Time (\mus)')
ylabel('Amplitude (a.u.)')

%% save data and compare
P_20mu = P_end;

figure
plot(T*1E6,P_0mu,T*1E6,P_100mu)
axis([[T(1) T(end)]*1E6 [-1 1]])
xlabel('Time (\mus)')
ylabel('Amplitude (a.u.)')

%% thickness inversion
% save data in single variable
Pinvert(:,1) = P_0mu;
Pinvert(:,2) = P_20mu;

% move to frequency domain
nfft = 2^nextpow2(length(T));
U = fft(Pinvert,nfft);
U = U(1:nfft/2+1,:)/nfft;
freq = Fn*linspace(0,1,nfft/2+1)';

% graphical check
figure
plot(freq*1E-6,abs(U))
xlabel('Frequency (MHz)')
ylabel('Amplitude (a.u.)')
axis([0 3 0 5E-3])

% experimental constants
rho_solid = 1200;
v_solid = c0;
rho_fluid = rho_solid;
v_fluid = c1;

% frequency band to integrate over
flow = find(freq>=.8E6,1);
fhigh = find(freq>=1.3E6,1);
freqband = freq(flow:fhigh);

% transmission coef
H = (-250:0.5:250)*1E-6;
alpha = 2*pi*freqband*H/v_fluid;  % freq * thickness
Zr = rho_fluid*v_fluid/(rho_solid*v_solid);
rff = (Zr-1)/(Zr+1);
Trans = ((1-rff^2)*exp(-1i*alpha))./(1-rff^2*exp(-2*1i*alpha));

% objective function
Fun = sum(abs((U(flow:fhigh,2)*ones(size(H))-(U(flow:fhigh,1)...
    *ones(size(H)).*Trans)).^2),1);

% locate min
[~, hmin] = min(Fun);
disp(H(hmin)*1E6)

% plot objective function
figure
plot(H*1E6,Fun,[H(hmin) H(hmin)]*1E6,[-1 1],':k')
axis([[H(1) H(end)]*1E6 0 3E-3])

%% old computation
% initial conditions
p = zeros(length(T),length(X));
q = zeros(length(T),length(X));

% boundary conditions
w = ricker(f0,dt);
p(1:length(w),1) = w;

% time n=2 with Lax-wendroff
for nn = 2:2
    for j = 2:length(X)-1
        p(nn,j) = p(nn-1,j) - C(j)/2 * (q(nn-1,j+1)-q(nn-1,j-1)) + C(j)^2 /2 *(p(nn-1,j-1)-2*p(nn-1,j)+p(nn-1,j+1)) ;
        q(nn,j) = q(nn-1,j) - C(j)/2 * (p(nn-1,j+1)-p(nn-1,j-1)) + C(j)^2 /2 *(q(nn-1,j-1)-2*q(nn-1,j)+q(nn-1,j+1)) ;
    end
end
figure

% time n>2 with Leap-Frog
for nn = 3:3 %length(t)
    for j = 2:length(X)-1
        p(nn,j) = 2*p(nn-1,j)-p(nn-2,j)+C(j)^2 * (p(nn-1,j+1)+p(nn-1,j-1)-2*p(nn-1,j));
    end
    
    % plotting at every kk step
    kk = 10;
    if mod(nn,kk) == 0
        plot(X,p(nn,:),'r',[x0 x0],[-1 1]*2,':k')
        title(['iteration n ' num2str(nn)])
        %axis([x(1) x(length(x)) -1.1 1.1])
        axis([X(1) x0*2 -1.1 1.1])
        xlabel('x')
        ylabel('p')
        pause(0.01)
    end
end