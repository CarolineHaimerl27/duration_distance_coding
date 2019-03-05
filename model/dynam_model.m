function y=dynam_model(t,x)

                        % N = x(1); tauR = x(4); 
Ipar=x(2); f_theta = x(3); I_theta=x(5);
global N Ufix tauR tauF W threshii


% J1 = 70; % strength of map-specific interaction (excitatory strength), given
% J0 = 50; % uniform feedback inhibition (inhibitory strength), given
% delta = 0.2; % strength of asymmetry in connections (in rad), given

%tauR = 0.1*tau_sc; % 0.1 time constant of recovery from synaptic depression, given
% tauF = 0.3; % 0.3 facilitation recovery time constant (in sec), given
tau = 10/1000; % tau is 10 millisecond, given
% external constant input
I = ones(1, N)*Ipar;
% input during break
% th = 1:N; theta = th/N*2*pi; % indizes of neural units
thresh = linspace(0,0.02,N)+threshii; % in Wang paper denoted as g, 
%thresh = [thresh(15:end) thresh(1:14)];
% Ufix = 0.15;
% synaptic weights
% abdiff = repmat(theta', 1, N)-repmat(theta, N, 1);
% W0 = exp(-(abdiff-delta).^2) +...
%     exp(-(abdiff-delta+2*pi).^2)  +...
%     exp(-(abdiff-delta-2*pi).^2);
% W = J1*W0' - J0; W(diag(ones(length(W),1))==1)=0;
% Noise
% Noise = NoisePwr*randn(N,1);

X = x(6:(N+5)); %
U = x(N+6:(end-N)); % if constant: Ufix; %
M = x(2*N+6:end);

IR = sum(W.*((M.*X.*U)*ones(1,N)))./(N-1); % recurrent contribution
IE = I+I_theta*sin(2*pi*f_theta*t);%+Noise'; % external input: constant and theta
x0 = (1-X)./tauR- U.*X.*M;
u0 = (Ufix-U)./tauF+Ufix.*(1-U).*M; % change in facilitation u
Iall = (IR+IE)-thresh;
Iall(Iall<=0)=0; %Iall(Iall>0)=15;
m0 = (-M+Iall')/tau; % change in firing rate

y=[zeros(5,1);x0; u0; m0];%x0; u0; m0];%  zeros(N,1); m0];%x0; zeros(N,1); m0];

end