function [t,x] = Basic_Model(parameters,x0,time_data)

global par W0 I Ipar superFlag


par = parameters;

N = par(1);
delta = par(end);
th = 1:N; theta = th/N*2*pi; % indizes of neural units

abdiff = repmat(theta', 1, N)-repmat(theta, N, 1);

W0 = exp(-(abdiff-delta).^2) +...
    exp(-(abdiff-delta+2*pi).^2)  +...
    exp(-(abdiff-delta-2*pi).^2);

opts = odeset('Events',@MyEvents);
[t,x,te,ye]=ode45(@dynam_model,time_data,x0,opts); % Runge Kutta integration

if(~isempty(te))
    x=0;
    t=0;
end

