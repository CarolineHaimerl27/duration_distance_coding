%%%%%%%% call firing rate model %%%%%%%
%%%%%%%%%%%%%%% + graphs %%%%%%%%%%%%%%

%% define parameters
global N Ufix tauR tauF W threshii
N       = 100;      % number of units
tend    = 3;        % length of trial in seconds
dt      = 1e-4;     % time step
Ufix    = 0.15;     % fixed point of facilitation U
tauR    = 0.3;      % time constant of depression
tauF    = 0.3;      % time constant of facilitation
Mscale  = 1;        % scaling the starting (baseline) firing rate
thetaon = 1;        % if =1, there is theta input
J1 = 70;            % strength of map-specific interaction (excitatory strength)
J0 = 50;            % uniform feedback inhibition (inhibitory strength)
delta = 0.2;        % strength of asymmetry in connections (in rad)
th = 1:N; theta = th/N*2*pi; % indizes of neural units
abdiff = repmat(theta', 1, N)-repmat(theta, N, 1);
tbreak  = 0;        % if >0 simulates an input break
threshii = 0;       % threshold for the nonlinearity
INPC    = 0;        % constant input
FREQ    = 9;        % theta input frequency
AMP     = 18;       % theta input amplitude, for slope over amplitude plot
                    % set AMP = 5:2:45; 
modelspec = 'dynam_model';
%% simulation
clear SIM;
count=1; 
if sum(INPC>0)>0, thetaon=0; end
 
% connectivity matrix
W0 = exp(-(abdiff-delta).^2) +...
    exp(-(abdiff-delta+2*pi).^2)  +...
    exp(-(abdiff-delta-2*pi).^2);
W = J1*W0' - J0; W(diag(ones(length(W),1))==1)=0;

% run simulation over different parameter ranges
for jj=1:length(AMP) 
    for ii=1:length(FREQ) 
        for cc=1:length(INPC)
            tic
            SIM(count) = FRModel(INPC(cc), tend, dt, Mscale, thetaon,... 
                                FREQ(ii), AMP(jj), tbreak, modelspec);
            toc
            count = count+1;
            drawnow 
            %close all
        end
    end
end

%% facilitation
ind = 1;
figure;
imagesc(SIM(ind).Udyn')
title('facilitation')
colorbar

%% depression
figure;
imagesc(SIM(ind).Xdyn')
title('depression')
colorbar

%% example neuron
figure
subplot(3,1,1)
hold
plot([0, max(SIM.t)], [Ufix, Ufix])
plot(SIM(ind).t, SIM(ind).Udyn(:,1))
title('facilitation')
subplot(3,1,2)
hold
plot([0, max(SIM.t)], [1, 1])
plot(SIM(ind).t, SIM(ind).Xdyn(:,1))
title('depression')
subplot(3,1,3)
hold
plot(SIM(ind).t, INPC+AMP(ind)*sin(2*pi*FREQ*SIM(ind).t))
plot(SIM(ind).t, SIM(ind).M(:,1))
title('firing rate')

%% slope over input amplitude

figure; 
hold on
plot(AMP,cat(1,SIM.slope))


