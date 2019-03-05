function [thslope, thshift] = ThSlope(t, f_theta, N, M, dt)

NCycles = round(max(t)*f_theta);
Fmax    = nan(NCycles,N);
th      = max(M(:))/10; %threshold on activity
ThStart = nan(1,NCycles);
ThLength = nan(1,NCycles);
thslope = nan(1,NCycles);
for i = 1:NCycles
    tmp     = 1+round(1/f_theta/dt*(i-1)):round(1/f_theta/dt*i);
    [a,b]   = max(M(tmp(tmp<size(M,1)),:));
    %Fmax: position of max firing rate in theta cycle (rad)
    Fmax(i, a>th) = b(a>th)*f_theta*dt*2*pi;
    %Theta sequence slope, length and shift
    ThStarttmp = find(diff(isnan(Fmax(i,:))) == -1,1);
    ThEndtmp    = find(diff(isnan(Fmax(i,:))) == 1,1);%First neuron of theta sequence
    if ~isempty(ThStarttmp) || ~isempty(ThEndtmp) %Is there a biginning or an end ?
        if ~isempty(ThStarttmp)
            ThStart(i) = ThStarttmp;
        else
            ThStart(i) = 1;
        end
        ThLength(i) = sum(~isnan(Fmax(i,:)));
        %Unwrap in neuron space
        tmp2        = circshift(Fmax(i,:),[0 -ThStart(i)]);
        %Fit theta sequence
        if sum(isnan(tmp2)==0)>3; % minimum of points to perform regression on
            X           = robustfit(tmp2,(1:N));
            thslope(i)  = X(2)*2*pi*f_theta; % (Neuron/s)
        end
    end
end
thshift = diff(unwrap(ThStart*2*pi/N))*N/2/pi;