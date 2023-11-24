% Here are some example scripts to simulate the AP model of:
% Aliev, Rubin R., and Alexander V. Panfilov. "A simple two-variable model of cardiac excitation." Chaos, Solitons & Fractals 7.3 (1996): 293-301.
%
% The equations here are based on the simplified form from:
% http://www.ibiblio.org/e-notes/html5/ap.html
%
% AEC 10 July 2015

%% First explore the model for a single cell.
t0 = 0;  % Integration window
tfinal = 55;
u0 = 0.08; % Initial conditions for the two state-variables: u, v.
v0 = 0;
a = 0.05;  % Two system parameters: a, k
k = 8;
c = 0.03;  % Driver current

% Simulation for a single oscillator
[t, y] = aliev(t0, tfinal, u0, v0, a, k, c);

figure(1); clf
plot(t, y)
legend('u (voltage)', 'v (hidden variable)')
xlabel('Time (A.U.)')

% explore how beat rate depends on drive current, c
clear t y
c = ((0:50)/50) * 0.05 + 0.01;
for j = 1:length(c);
    [t{j}  y{j}] = aliev(t0, tfinal, u0, v0, a, k, c(j));
end;
    
figure(1)
for j = 1:length(c);
    plot(t{j}, y{j}(:,1) + j/10);
    hold all
end;
hold off;
% beat rate increases with steady current between c = 0.025 and 0.039.

% explore how beat rate depends on parameter 'a'
clear t y
c = 0.02;
a = ((0:50)/50) * 0.07;
for j = 1:length(a);
    [t{j}  y{j}] = aliev(t0, tfinal, u0, v0, a(j), k, c);
end;
    
figure(1)
for j = 1:length(a);
    plot(t{j}, y{j}(:,1) + j/10);
    hold all
end;
hold off;

% explore how beat rate depends on parameter 'k'
clear t y
c = 0.02;
a = 0.04;
k = ((0:50)/50)*10;
for j = 1:length(k);
    [t{j}  y{j}] = aliev(t0, tfinal, u0, v0, a, k(j), c);
end;
    
figure(1)
for j = 1:length(k);
    plot(t{j}, y{j}(:,1) + j/10);
    hold all
end;
hold off;

% Make a summary plot
c = [0 .02 .04 .06];
a = [0 .02 .04 .06];
k = ((0:50)/50)*10;
for p = 1:4;
    figure(1); clf
    for q = 1:4;
        clear t y
        for j = 1:length(k);
            [t{j}  y{j}] = aliev(t0, tfinal, u0, v0, a(q), k(j), c(p));
        end;
        subplot(2,2,q)
        for j = 1:length(k);
            plot(t{j}, y{j}(:,1) + j/10);
            hold all
        end;
        title(['a = ' num2str(a(q))])
        text(10, 6, ['k = 0:0.2:10; c = ' num2str(c(p))])
        ylim([0 7])
        hold off;
    end;
    saveas(gca, ['Aliev model c = ' num2str(c(p)) '.fig']);
    saveas(gca, ['Aliev model c = ' num2str(c(p)) '.png']);
end;


%%
% Now work with a 2D model.

% Simulate islands of different sizes to look at the effect of island size
% on spiking patterns.  Give each cell heterogeneous drive currents
nx = [2 3 5 10];  % nx = [2 5 10 20] % for the first attempt
ny = [2 3 5 10];  % ny = [2 5 10 20] % for the first attempt
phi = [.05 .2];   % phi = [0 .1 .2 .4 .8] % for the first attempt
t0 = 0;
tfinal = 60;

tic
for LIdx = 1:length(nx);
    u0 = 0.08*zeros(ny(LIdx), nx(LIdx));
    v0 = zeros(ny(LIdx), nx(LIdx));
    a = 0.05*ones(ny(LIdx), nx(LIdx));
    k = 8*ones(ny(LIdx), nx(LIdx));
    for phiIdx = 1:length(phi);
        for j = 1:50;
            c = 0.032*ones(ny(LIdx), nx(LIdx)) + 0.004*randn(ny(LIdx), nx(LIdx));  % drive current
            'j phi L'
            [j phi(phiIdx) nx(LIdx)]
            [t, uDat, vDat] = aliev2D(t0, tfinal, u0, v0, a, k, c, nx(LIdx), ny(LIdx), phi(phiIdx));
            intens{j,phiIdx,LIdx} = squeeze(mean(mean(uDat)));
            tau{j,phiIdx,LIdx} = t;
        end;
    end;
end;
toc
save island_size_sim2.mat intens tau nx ny phi t0 tfinal u0 v0 a k c

load('island_size_sim2.mat')


% It's hard to work with nonuniform time axes, so make a uniform time axis:
tAx = 0:.04:60;
intensI = zeros(length(tAx), 50, length(phi), length(nx));

for LIdx = 1:length(nx);
    for phiIdx = 1:length(phi);
        for j = 1:50;
            intensI(:,j,phiIdx,LIdx) = interp1(tau{j,phiIdx,LIdx}, intens{j,phiIdx,LIdx}, tAx);
        end;
    end;
end;

% Look at a sampling of the results
for LIdx = 1:length(nx);
    for phiIdx = 1:length(phi);
        subplot(length(nx),length(phi),(LIdx - 1)*length(phi) + phiIdx)
        for j = 1:2:50;
            plot(tAx, j/50+intensI(:,j, phiIdx, LIdx)); hold all
            title(['L = ' num2str(nx(LIdx)) '; phi = ' num2str(phi(phiIdx))])
%           plot(tau{j,phiIdx, LIdx}); hold all
        end;
        hold off;
    end;
end;
saveas(gca, 'Example traces.fig')
saveas(gca, 'Example traces.png')

% Use the time interval between the last and the second-to-last peak as a
% measure of the period
peaksDat = zeros(50,length(phi),length(nx));
for LIdx = 1:length(nx);
    for phiIdx = 1:length(phi);
        for j = 1:50;
            [~,tmp] = findpeaks(intensI(:,j,phiIdx,LIdx));
            peaksDat(j,phiIdx,LIdx) = tAx(tmp(end) - tmp(end-1));            
        end;
    end;
end;

% Summary statistics: mean and standard deviation in beat rate vs island
% size
figure(2); clf
subplot(1,2,1)
meandat = squeeze(mean(peaksDat(:,:,:), 1))';
stddat = squeeze(std(peaksDat(:,:,:),[], 1))';
for phiIdx = 1:length(phi);
    errorbar(nx, meandat(:,phiIdx), stddat(:,phiIdx)/sqrt(50))
    hold all
end;
hold off
xlabel('Island edge length')
ylabel('Mean beat period')
ylim([15 17])
subplot(1,2,2)
plot(nx, stddat)
xlabel('Island edge length')
ylabel('Standard deviation beat period')
saveas(gca, 'mean and std vs island size.fig')
saveas(gca, 'mean and std vs island size.png')




%% Other games with the 2D system
% Try to get a spiral wave
nx = 20;
ny = 20;
phi = 0.2;
c = 0.01*ones(ny, nx);  % drive current
% Impose a spiral starting wave
theta = 0:.01:(pi);
xcoord = cos(2*theta).*theta/(pi);
ycoord = sin(2*theta).*theta/(pi);
xcoord = round(nx/2 + nx/2.5*xcoord);
ycoord = round(ny/2 + ny/2.5*ycoord);
for j = 1:length(theta); c(ycoord(j), xcoord(j)) = 0.06; end;
figure(12)
imshow(c, [], 'InitialMagnification', 'fit')

t0 = 0;
tfinal = 10;
u0 = zeros(ny, nx);
dtheta = -.4;
theta = (0:.01:(pi));
xcoord = cos(2*(theta + dtheta)).*theta/(pi);
ycoord = sin(2*(theta + dtheta)).*theta/(pi);
xcoord = round(nx/2 + nx/2.5*xcoord);
ycoord = round(ny/2 + ny/2.5*ycoord);
for j = 1:length(theta); u0(ycoord(j), xcoord(j)) = 0.08; end;
figure(13)
imshow([c u0], [], 'InitialMagnification', 'fit')

v0 = zeros(ny, nx);
a = 0.02;
k = 10;

% Break the simulation into two parts, first with the external current
% imposed
[t1, uDat1, vDat1] = aliev2D(t0, tfinal, u0, v0, a, k, c, nx, ny, phi);

% then with a uniform conditions
c = 0.02*ones(ny, nx);
tfinal = 60;
[t2, uDat2, vDat2] = aliev2D(t0, tfinal, uDat1(:,:,end), vDat1(:,:,end), a, k, c, nx, ny, phi);
t2 = t2 + max(t1);
uDat = cat(3, uDat1, uDat2);
t = [t1; t2];
nT = length(t)

figure(3); clf
for j = 1:nT;
    imshow(uDat(:,:,j), [0, 1], 'InitialMagnification', 'fit');
    title(num2str(j))
    pause(.01);
end;


%% Play with point sources in 1-D
nx = 20;
ny = 1;
phi = 0.2;
c = 0.02*ones(ny, nx);  % drive current


% Impose a spiral starting wave
c(1,10) = 0.05;

t0 = 0;
tfinal = 60;
u0 = zeros(ny, nx);
v0 = zeros(ny, nx);
a = 0.05;
k = 8;

[t, uDat, vDat] = aliev2D(t0, tfinal, u0, v0, a, k, c, nx, ny, phi);
% c = 0.00*ones(ny, nx);
% tfinal = 60;
% [t2, uDat2, vDat2] = aliev2D(t0, tfinal, uDat1(:,:,end), vDat1(:,:,end), a, k, c, nx, ny, phi);
% t2 = t2 + max(t1);
% uDat = cat(3, uDat1, uDat2);
% t = [t1; t2];
nT = length(t)
uDat = squeeze(uDat)';
tAx = ones(20,1)*(1:1000)*60/1000;
Xmat = ones(nT,1)*(1:20);
Tmat = t*ones(1,20);
uInterp = interp2(Tmat', Xmat', uDat', (1:.1:60)', (1:20))';
figure(3); clf
subplot(2,1,1);
pcolor(uInterp); shading 'interp'
subplot(2,1,2);
plot(uInterp(:,8))


figure(3); clf
for j = 1:nT;
    imshow(uDat(:,:,j), [0, 1], 'InitialMagnification', 'fit');
    title(num2str(j))
    pause(.01);
end;

%%
%% Other games with the 2D system
% Are there any conditions where a single cell can launch a wave?
nx = 20;
ny = 20;
phi = 0.2;
c = 0.0*ones(ny, nx);  % drive current

t0 = 0;
tfinal = 60;
u0 = zeros(ny, nx);
u0(10,10) = 0.4;

v0 = zeros(ny, nx);
a = 0.04;
k = 8;

[t, uDat, vDat] = aliev2D(t0, tfinal, u0, v0, a, k, c, nx, ny, phi);
nT = length(t)

figure(3); clf
for j = 1:nT;
    imshow(uDat(:,:,j), [0, 1], 'InitialMagnification', 'fit');
    title(num2str(j))
    pause(.01);
end;

%% 2D spiral wave
% added by H McNamara, 16 July 2015

nx = 512;
ny = 512;
phi = 0.05; % if coupling too high, cannot block wave to get spiral pattern
c = 0.0*ones(ny, nx);  % drive current

t0 = 0;
tfinal = 300;
u0 = zeros(ny, nx);
u0(256,257:end) = 1;

% create boundary for directional propagation
v0 = zeros(ny, nx);
v0(257,257:end) = 1;
a = 0.04;
k = 8;

tic
[t, uDat, vDat] = aliev2D(t0, tfinal, u0, v0, a, k, c, nx, ny, phi);
nT = length(t);
toc

% check final frame for spiral pattern
imshow(uDat(:,:,nT),[0,1]);

% watch movie
figure(3); clf
for j = 1:50:nT;
    imshow(uDat(:,:,j), [0, 1], 'InitialMagnification', 'fit');
    title(num2str(j))
    pause(.01);
end;
