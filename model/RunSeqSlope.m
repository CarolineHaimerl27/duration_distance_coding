function slope = RunSeqSlope(M, N, t, optth, f_theta)
mfname = mfilename;
if nargin < 3 || isempty(f_theta)
    warning('%s: no theta freq input, estimation point by point', mfname);
    optth=0;
end

% find slope (careful, no knockout-criteria for slope here!!!)
if optth ==0;
% use only 80% of simulation - discard beginning
M = M(round(size(M,1)*0.2):end,:);
[Mmaxtmp,Ind] = max(M,[],2);      % Maxima of firing rate at each point in time
Ind(Mmaxtmp < max(M(:))/10) = NaN;              % Remove points that are lower than 10% of maximum
Ind = unwrap(Ind/N*2*pi); %Unwrap due to circular topology (conversion in rad)
Ttim = round(size(M,1)*0.2);
ttmp = t(1:Ttim); Indtmp = Ind((end-Ttim+1):end); %Select stable solution
X = polyfit(ttmp(~isnan(Indtmp)),Indtmp(~isnan(Indtmp)),1);
slope = X(1)*N/2/pi; % Convert back to neurons (Slope = neurons/s)
figure; imagesc(t, 1:N, M')
    hold on
    plot(t, t*slope, 'r')
    
else
%% slope detection if there is theta
theta = sin(2*pi*f_theta*t);
% find theta around peaks --> maximum input --> maximum theta blob
[~, ptheta] = findpeaks(1.*(theta.*(abs(theta)<0.1)), 'MinPeakHeight', 0);
% group time points into groups of theta blobs
[~, b] = histc(t, t(ptheta));
% determine the sequence of maximum active neurons in each theta blob
mfr = []; tmfr = []; mfrm = []; tmfrm = []; 
figure; imagesc(t, 1:100, M'); 
hold on
if max(b) > 3; % if there are more than 3 theta cycles
    for jj=3:2:max(b)
        [ax, x] = max(M(b==jj,:)');
        mfr = [mfr x];
        tmfr = [tmfr t(b==jj)'];
        mfrm = [mfrm x(ax==max(ax))];
        tmfrm = [tmfrm median(t(b==jj))];
         plot(t(b==jj),x, 'b.')
         plot(median(t(b==jj)),x(ax==max(ax)), 'ro')
    end
    % exclude the minimum and the maximum value since in the beginning
    % there might be a fake-slope
    mfrm2 = mfrm; mfrm2(mfrm==max(mfrm)|mfrm==min(mfrm))=[];
    % test if slope exists by checking if peaking neurons span the whole
    % neuron space --> check range in peaking neurons, should be at least
    % 70%
    if range(mfrm2)>(N*0.70);
                            %     the distribution of peaking
                            %     % neurons is somewhat uniformly or peaking
                            % hn = histc(mfrm, 1:round(N*0.3):(N+N*0.3));
                            %     if range(hn)> max(hn)*0.5;
                            %         warning('no slope detected')
                            %         slope = NaN;
        % unwrap due to circular topology
        MFR = unwrap(mfr/N*2*pi); %Unwrap due to circular topology (conversion in rad)
        MFRm = unwrap(mfrm/N*2*pi); %Unwrap due to circular topology (conversion in rad)
        % if mfrm does not span a sufficient part of the cell space - no
        % slope can be defined! unwrap only takes care of jumps>pi, but for
        % jumps that are <pi can also be problematic!
        % set a boundary at jumps of 10% of 2pi
        if median(abs(diff(MFRm)))>(2*pi*0.35)
            warning('jumps between max neuron in theta cycles larger than 20% of N --> no slope detected')
            slope = NaN;
        else
            C = polyfit(tmfr,MFR, 1);
            C2 = polyfit(tmfrm, MFRm, 1); 
            figure; plot(tmfrm, MFRm*N/2/pi, 'rx'); hold on; 
                    plot(tmfr, MFR*N/2/pi, 'bx'); 
                    plot(tmfr, polyval(C*N/2/pi, tmfr), 'b');
                    plot(tmfr, polyval(C2*N/2/pi, tmfr), 'r');
                    legend('mean', 'all')
             if abs(C2(1)-C(1))>(max(abs(C2(1)), abs(C(1)))*0.1);
                warning('slope detection with max neuron vs max neuronS in theta cycle varies by more than 10% of slope, use max neuron method')
%                 disp(C(1)*N/2/pi)
%                 disp(slope)
            end
            % how does the fit perform? use r^2
            ei = polyval(C2, tmfrm)-MFRm;
            ym = mean(MFRm);
            SST = sum((MFRm-ym).^2);
            SSR = sum(ei.^2);
            R2 = 1-SSR/SST;
            if R2<0.95;
                 warning('Fit criteria r^2<0.95, sequence slope not linear!')
                 slope = NaN;
            else
                % Convert back to neurons (Slope = neurons/s)
                slope = C2(1)*N/2/pi; 
            end
                 
                 
                 
        end
    else
        warning('max neuron in theta cycles concentrated around an area<70% of N, no slope detected')
        slope = NaN;
    end
else
        warning('%s: cannot estimate slope, not enough theta cycles', mfname);
        slope = NaN;
end

end