function y=dynam_model(t,x)
     
global par W0 superFlag fTheta

N=par(1);
J1=par(2);
J0=par(3); 
tauF=par(4); 
tauR = par(5);
tau=par(6);
Ufix=par(7); 
I = par(8);
fTheta = par(9);


if(max(x(end/2+1:end))>200)  
    superFlag = 1;
end

W =  J1*W0' - J0;

M = x(1:N);
U = x(N+1:2*N);
D = x(2*N+1:end);

IR = sum(W.*((M.*U.*D)*ones(1,N)))./(N-1); % recurrent contribution

% These two lines are for simulating constant or theta modulated input
%IE = I;
IE = I*sin(2*pi*fTheta*t); % external input: constant and theta

u0 = (Ufix-U)./tauF+Ufix.*(1-U).*M;
d0 = (1-D)./tauR-U.*D.*M;

Iall = (IR+IE);        
Iall(Iall<=0)=0;
               
m0 = (-M+Iall')/tau;
y=[m0;u0;d0];        

end