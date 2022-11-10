close all;

% Constants
EK = -82; ENA = 45; ELEAK = -59.4001; %milli volts
VREST = -70; %milli volts
GKBAR = 36; GNABAR = 120; GLEAK = 0.3; %milli siemens /cm^2
CMEM = 1; %micro farad
TIMESTEP = 0.001; %as our time unit is in milliseconds, this corresponds to 1 microsecond

SIMULATION_TIME = 50; %milli seconds
t = 0:TIMESTEP:SIMULATION_TIME;

% Initializations
GACH = 0; %initially closed
Vmem = zeros(1,length(t));
m = zeros(1,length(t));
n = zeros(1,length(t));
h = zeros(1,length(t));
gk = zeros(1,length(t));
gna = zeros(1,length(t));

Vmem(1) = VREST; %initially the cell is at rest
V = Vmem(1) - VREST; 
alpha_m = 0.1*(25-V)/(exp(0.1*(25-V))-1);
alpha_n = 0.01*(10-V)/(exp(0.1*(10-V))-1);
alpha_h = 0.07*exp(-V/20);
beta_m = 4*exp(-V/18);
beta_n = 0.125*exp(-V/80);
beta_h = 1/(exp(0.1*(30-V))+1);

% I assumed m_inf etc. corresponds to steady-state values
% so picked those as the initial values as the cell is at rest initially
m(1) = alpha_m/(alpha_m + beta_m);
n(1) = alpha_n/(alpha_n + beta_n);
h(1) = alpha_h/(alpha_h + beta_h);

for i = 1:length(t)-1
    V = Vmem(i) - VREST;
    alpha_m = 0.1*(25-V)/(exp(0.1*(25-V))-1);
    alpha_n = 0.01*(10-V)/(exp(0.1*(10-V))-1);
    alpha_h = 0.07*exp(-V/20);
    beta_m = 4*exp(-V/18);
    beta_n = 0.125*exp(-V/80);
    beta_h = 1/(exp(0.1*(30-V))+1);

    gk(i) = GKBAR * (n(i)^4);
    gna(i) = GNABAR * (m(i)^3) * h(i);

    ina = (gna(i) + GACH) * (Vmem(i) - ENA);
    ik = (gk(i) + GACH) * (Vmem(i) - EK);
    il = GLEAK * (Vmem(i) - ELEAK);
    iion = -(ina + ik + il);

    dm = alpha_m * (1-m(i)) - beta_m * m(i);
    dn = alpha_n * (1-n(i)) - beta_n * n(i);
    dh = alpha_h * (1-h(i)) - beta_h * h(i);

    Vmem(i+1) = Vmem(i) + TIMESTEP*iion/CMEM;
    m(i+1) = m(i) + TIMESTEP*dm;
    n(i+1) = n(i) + TIMESTEP*dn;
    h(i+1) = h(i) + TIMESTEP*dh;

    if i > 1/TIMESTEP
        GACH = 0.89; % set to 0 for no stimulation / 0.89-1ms / 0.24-1.5ms / 0.1392-2ms / 0.09-3ms
    end
    if i > 2/TIMESTEP
        GACH = 0;
    end
    if i > 12/TIMESTEP % the second pulse starts to appear after 2ms of the first but is not fully formed until 5.3 ms passes
        GACH = 0.89;
    end
    if i > 13/TIMESTEP
        GACH = 0;
    end
end
plot(t,Vmem)
hold on
plot(t,VREST*ones(1,length(t)),'r--')
legend('V_{mem}', 'V_{rest}')
ylabel('Voltage (mV)')
xlabel('time (ms)')
title('Potential across the Cell Membrane for Simulated Neuron ')

figure
plot(t,gk,'b')
hold on
plot(t,gna,'r')
hold on
legend('g_{K}', 'g_{Na}')
ylabel('Conductance (milli siemens/cm^{2})')
xlabel('time (ms)')
title('g_{K} and g_{Na} for Simulated Neuron')

conductances = [0.89 0.24 0.1392 0.09];
durations = [1 1.5 2 3];
figure
plot(durations,conductances)
ylabel('Conductances (milli siemens/cm^{2})')
xlabel('Stimuli durations (ms)')
