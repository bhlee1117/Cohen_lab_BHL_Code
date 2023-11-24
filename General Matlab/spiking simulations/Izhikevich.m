% Simulation of the Izhikevich model

clear
nSteps = 10000;
dt = 0.1; % ms
a = 0.02; %
b = 0.2; %
c = -65; % 
d = 4; % 

v = zeros(nSteps, 1);
v(1) = -65; % Initial values of v
u = zeros(nSteps, 1);
u(1) = b*v(1); % Initial values of u
firings=[]; % spike timings
I = zeros(nSteps, 1);
nStart = 101;
% Linear current ramp
I(nStart:end) = 50*((nStart:nSteps)-nStart)/(nSteps - nStart);

for j=1:nSteps - 1
    v(j+1) = v(j) + dt*(0.04*v(j)^2 + 5*v(j) + 140 - u(j) + I(j)); % 
    u(j+1) = u(j) + dt*a*(b*v(j+1) - u(j)); 
    if v(j + 1) >= 30;
        firings = [firings j];
        v(j) = 30;
        v(j+1) = c;
        u(j+1) = u(j+1)+d;
    end;
end;

t = (1:nSteps)*dt;
plot(t, v, t, u)

% % Look at the relation between v and dv/dt
% v1 = -100:100;
% vdot = 0.04*v1.^2 + 5*v1 + 140;
% plot(v1, vdot)
% xlabel('Voltage')
% ylabel('dV/dt')
%%

% Try random values of the parameters:
% Under current clamp drive
v = zeros(nSteps, 1);
v(1) = -65; % Initial values of v
u = zeros(nSteps, 1);
u(1) = b*v(1); % Initial values of u
firings=[]; % spike timings

% Under optogenetic drive:
v2 = zeros(nSteps, 1);
v2(1) = -65; % Initial values of v
u2 = zeros(nSteps, 1);
u2(1) = b*v2(1); % Initial values of u
firings2=[]; % spike timings

t = (1:nSteps)*dt;
while (1)
    I = 5 + 5*randn(nSteps, 1);
    alpha = I/-65; % photoconductance.  When v = -65 mV, the photoconductance and current are matched.
    a = .025 + .05*rand;
    b = .1 + .2*rand;
    c = -65 + 10*rand;
    d = 0 + 6*rand;
    u(1) = b*v(1);
    u2(1) = b*v2(1);
    % Current injection
    for j=1:nSteps - 1
        v(j+1) = v(j) + dt*(0.04*v(j)^2 + 5*v(j) + 140 - u(j) + I(j)); % 
        u(j+1) = u(j) + dt*a*(b*v(j+1) - u(j)); 
        if v(j + 1) >= 30;
            firings = [firings j];
            v(j) = 30;
            v(j+1) = c;
            u(j+1) = u(j+1)+d;
        end;
    end;
    
    % ChR stim
    for j=1:nSteps - 1
        v2(j+1) = v2(j) + dt*(0.04*v2(j)^2 + 5*v2(j) + alpha(j)*v2(j) + 140 - u2(j)); %   conductance change!
        u2(j+1) = u2(j) + dt*a*(b*v2(j+1) - u2(j)); 
        if v2(j + 1) >= 30;
            firings2 = [firings2 j];
            v2(j) = 30;
            v2(j+1) = c;
            u2(j+1) = u2(j+1)+d;
        end;
    end;
    subplot(2,1,1)
    plot(t, v)
    xlabel('Time (ms)')
    ylabel('Voltage (mV)')
    title(['a = ' num2str(a) ' b = ' num2str(b) ' c = ' num2str(c) ' d = ' num2str(d)]);
    subplot(2,1,2)
    plot(t, v2)
    xlabel('Time (ms)')
    ylabel('Voltage (mV)')
    pause
end;

% Compare # spikes under ChR drive vs. I-clamp drive
nRuns = 400;
nSpikes = zeros(nRuns, 2);
for k = 1:nRuns;
    I = 5 + 5*randn(nSteps, 1);
    alpha = I/-65; % photoconductance.  When v = -65 mV, the photoconductance and current are matched.
    a = .025 + .05*rand;
    b = .1 + .2*rand;
    c = -65 + 10*rand;
    d = 0 + 6*rand;
    u(1) = b*v(1);
    u2(1) = b*v2(1);
    % Current injection
    firings = [];
    for j=1:nSteps - 1
        v(j+1) = v(j) + dt*(0.04*v(j)^2 + 5*v(j) + 140 - u(j) + I(j)); % 
        u(j+1) = u(j) + dt*a*(b*v(j+1) - u(j)); 
        if v(j + 1) >= 30;
            firings = [firings j];
            v(j) = 30;
            v(j+1) = c;
            u(j+1) = u(j+1)+d;
        end;
    end;
    nSpikes(k,1) = length(firings);
    
    % ChR stim
    firings2 = [];
    for j=1:nSteps - 1
        v2(j+1) = v2(j) + dt*(0.04*v2(j)^2 + 5*v2(j) + alpha(j)*v2(j) + 140 - u2(j)); %   conductance change!
        u2(j+1) = u2(j) + dt*a*(b*v2(j+1) - u2(j)); 
        if v2(j + 1) >= 30;
            firings2 = [firings2 j];
            v2(j) = 30;
            v2(j+1) = c;
            u2(j+1) = u2(j+1)+d;
        end;
    end;
    nSpikes(k,2) = length(firings2);
    k
end;
figure(2); clf
plot(nSpikes(:,1), nSpikes(:,2), 'x')
xlabel('I_{clamp}')
ylabel('ChR')
