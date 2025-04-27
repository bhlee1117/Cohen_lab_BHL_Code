function Vm=simulate_PSP(E_spk, I_spk, O_stim, V_rest, dt, T, E_opto, noise_std)
% Simulates full Vm(t) trace with passive membrane model and conductance-based inputs.
%
% Inputs:
%   E_spk, I_spk, O_stim - binary vectors (1 x T)
%   V_rest  - resting membrane potential (mV)
%   dt      - time step (s)
%   T       - total time (s)
%   E_opto  - reversal potential for optogenetic input (mV)

% --- Parameters ---
Cm = 200e-12;   % capacitance (Farads)
gL = 10e-9;     % leak conductance (Siemens)
EL = V_rest;    % leak reversal (mV)

% Time vector
time = 0:dt:T-dt;

% Reversal potentials (mV)
E_exc = 0;
E_inh = -70;

% Synaptic kinetics (tau in sec)
tau_rise_exc = 3.9e-3; tau_decay_exc = 10e-3;
tau_rise_inh = 2e-3; tau_decay_inh = 20e-3;
tau_rise_opto = 5e-3; tau_decay_opto = 15e-3;

% Conductance kernels (unitless, will be scaled)
t_kern = 0:dt:0.05;
gE_kernel = (exp(-t_kern/tau_decay_exc) - exp(-t_kern/tau_rise_exc));
gI_kernel = (exp(-t_kern/tau_decay_inh) - exp(-t_kern/tau_rise_inh));
gO_kernel = (exp(-t_kern/tau_decay_opto) - exp(-t_kern/tau_rise_opto));

% Normalize and scale (nS)
gE_max = 5e-9; gI_max = 10e-9; gO_max = 5e-10;
gE = conv(E_spk, gE_kernel, 'same') * gE_max;
gI = conv(I_spk, gI_kernel, 'same') * gI_max;
gO = conv(O_stim, gO_kernel, 'same') * gO_max;

% --- Simulate Vm ---
Vm = zeros(1, length(time));
Vm(1) = V_rest;

for t = 2:length(time)
    %g_total = gL + gE(t) + gI(t) + gO(t);
    I_total = ...
        gL * (EL - Vm(t-1)) + ...
        gE(t) * (E_exc - Vm(t-1)) + ...
        gI(t) * (E_inh - Vm(t-1)) + ...
        gO(t) * (E_opto - Vm(t-1));

    dVm = (I_total / Cm) * dt;
    noise = randn * noise_std;
    Vm(t) = Vm(t-1) + dVm + noise;
end

% --- Plot ---
clf;
subplot(4,1,1); plot(time, gE*1e9, 'r'); ylabel('g_E (nS)');
title(sprintf('EPSC: \\tau_{rise}=%.1f ms, \\tau_{decay}=%.1f ms, g_{max}=%.1f nS', ...
    tau_rise_exc*1e3, tau_decay_exc*1e3, gE_max*1e9));

subplot(4,1,2); plot(time, gI*1e9, 'b'); ylabel('g_I (nS)');
title(sprintf('IPSC: \\tau_{rise}=%.1f ms, \\tau_{decay}=%.1f ms, g_{max}=%.1f nS', ...
    tau_rise_inh*1000, tau_decay_inh*1000, gI_max*1e9));

subplot(4,1,3); plot(time, gO*1e9, 'color',[0 0.6 1]); ylabel('g_O (nS)');
title(sprintf('Opto: \\tau_{rise}=%.1f ms, \\tau_{decay}=%.1f ms, g_{max}=%.1f nS', ...
    tau_rise_opto*1e3, tau_decay_opto*1e3, gO_max*1e9));

subplot(4,1,4); plot(time, Vm, 'k'); ylabel('V_m (mV)'); xlabel('Time (s)');
title('Membrane Voltage');

sgtitle(sprintf('Vm Simulation | V_{rest}=%.1f mV, E_{opto}=%.1f mV, Noise=%.2f mV', ...
    V_rest, E_opto, noise_std));
end
