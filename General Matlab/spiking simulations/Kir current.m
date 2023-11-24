%% Inward Rectifier Current
clear
Ek = -80;
V = -80:.1:20;
GK1 = (.1:.02:1)';  % vary the strength of the Kir current
nG = length(GK1);  % # of GK1 values.

% model for the Kir current
a_k1=0.1./(1+exp(0.06*(V-Ek-200)));
b_k1=(3*exp(0.0002*(V-Ek+100))+exp(0.1*(V-Ek-10)))./(1+exp(-0.5*(V-Ek)));
xK1_inf=(a_k1./(a_k1+b_k1));
IK1 = GK1*(xK1_inf.*(V-Ek));

% Assume a purely ohmic leak
ILeak = V*.001;
Itot = IK1 + repmat(ILeak, [nG, 1]);
plot(V, Itot)

% Find the zero-crossings of Itot
Isign = sign(Itot);
idx1 = []; idx2 = [];
for j = 1:nG;
    idx1(j) = find(Isign(j,:) == 1, 1);
    idx2(j) = find(Isign(j,:) == -1, 1, 'last');
end;

% show the zero-crossings on a plot to check for accuracy
plot(V, Itot); hold all;
for j = 1:nG;
    plot([V(idx1(j)), V(idx2(j))], [Itot(j,idx1(j)), Itot(j,idx2(j))], 'r*');
end;
hold off;
xlabel('V (mV)')
ylabel('i (nA)')

% plot(1:nG, idx1, 1:nG, idx2)
idx1(idx1 > idx2) = idx2(idx1 > idx2);

% Look at the equilibrium voltages
plot(GK1, V(idx1), GK1, V(idx2))
xlabel('GK1')
ylabel('Equilibrium voltages (mV)')

% The state that wins has the bigger area under the curve.
% Look at the cumulative sum of the current
for j = 1:nG;
    bal(j) = sum(Itot(j,idx1(j):idx2(j)),2);
end;
bal(idx1 == idx2) = 0;
plot(GK1, bal)

plot(GK1, V(idx2), GK1, V(idx1)); hold all
plot([GK1(bal < 0); GK1(bal > 0)], [V(idx2(bal < 0)) V(idx1(bal > 0))], 'k--');
plot(GK1(bal == 0), V(idx1(bal == 0)), 'k*');
hold off
xlabel('GK1')
ylabel('Equilibrium voltages (mV)')
legend('Stable point 1', 'Stable point 2', 'Syncytium')
saveas(gca, 'Equilibrium voltage with Kir.fig')
saveas(gca, 'Equilibrium voltage with Kir.png')



%% Repeat, varying Ileak for constant GK1
clear
Ek = -80;
V = -80:.1:20;
GK1 = .5;  % strength of the Kir current

% model for the Kir current
a_k1=0.1./(1+exp(0.06*(V-Ek-200)));
b_k1=(3*exp(0.0002*(V-Ek+100))+exp(0.1*(V-Ek-10)))./(1+exp(-0.5*(V-Ek)));
xK1_inf=(a_k1./(a_k1+b_k1));
IK1 = GK1*(xK1_inf.*(V-Ek));

% Assume a purely ohmic leak
GLeak = (.0001:.0001:.005);
ILeak = GLeak'*V;
nL = size(ILeak, 1);
Itot = repmat(IK1, [nL, 1]) + ILeak;
plot(V, Itot)

% Find the zero-crossings of Itot
Isign = sign(Itot);
for j = 1:nL;
    idx1(j) = find(Isign(j,:) == 1, 1);
    idx2(j) = find(Isign(j,:) == -1, 1, 'last');
end;

% show the zero-crossings on a plot to check for accuracy
plot(V, Itot); hold all;
for j = 1:nL;
    plot([V(idx1(j)), V(idx2(j))], [Itot(j,idx1(j)), Itot(j,idx2(j))], 'r*');
end;
hold off;
xlabel('V (mV)')
ylabel('i (nA)')

plot(1:nL, idx1, 1:nL, idx2)
idx1(idx1 > idx2) = idx2(idx1 > idx2);  % fix the corner case

% Look at the equilibrium voltages
plot(GLeak, V(idx1), GLeak, V(idx2))
xlabel('GLeak')
ylabel('Equilibrium voltages (mV)')

% The state that wins has the bigger area under the curve.
% Look at the cumulative sum of the current
bal = [];
for j = 1:nL;
    bal(j) = sum(Itot(j,idx1(j):idx2(j)),2);
end;
bal(idx1 == idx2) = 0;
plot(GLeak, bal)


plot(GLeak, V(idx2), GLeak, V(idx1)); hold all
plot([GLeak(bal > 0) GLeak(bal < 0)], [V(idx1(bal > 0)) V(idx2(bal < 0))], 'k--');
plot(GLeak(bal == 0), V(idx1(bal == 0)), 'k*');
hold off
xlabel('GLeak')
ylabel('Equilibrium voltages (mV)')
legend('Stable point 1', 'Stable point 2', 'Syncytium')
saveas(gca, 'Equilibrium voltage with varying GLeak.fig')
saveas(gca, 'Equilibrium voltage with varying GLeak.png')



