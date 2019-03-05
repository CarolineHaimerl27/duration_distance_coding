function [ SIMc, Cs, cyc, dirdet ] = SLOPE( M, t, f_theta, I_theta, N, ws, opt)
% M is the firing rate of each neuron
% t gives the time steps (can vary if RK is used)
% f_theta is frequency
% I_theta is amplitude
% N is number of neurons
% ws is window size to determine the slope (not to small otherwise slope
%       fitting won't work
% opt is 1 if graphs are wanted

% find the most active neuron
T=length(M); M=M';
if I_theta>0;
signal = I_theta*sin(2*pi.*f_theta.*t);
                        % s2 = signal(T*0.05:end); 
 [~, cyc]=findpeaks(signal); cyc=[cyc' size(M,2)]; % determine theta cycles over peaks in oscillation
else
    cyc = 2:2:size(M,2);
end
 for i=1:length(cyc)-1; % get the neuron with the highest mean firing rate in this cycle
     [~,maxcell(i)]= max(mean(M(:,cyc(i):cyc(i+1)),2)); 
     maxval(i)= max(max(M(:,cyc(i):cyc(i+1)),[],2)); % get the maximum firing of the best cell
     maxtim(i)= round(median(cyc(i):cyc(i+1)));
     if length(f_theta)>1; mfreq(i) = abs(1/diff([t(cyc(i)),t(cyc(i+1))])); %median(f_theta(cyc(i)):f_theta(cyc(i+1))); 
     else mfreq(i) = f_theta; end
 end % get maximum mean firing cell per theta cycle
% figure; plot(t(maxtim), maxcell, 'o')

% determine direction of sequence (new!!!)
    if sum(diff(maxcell)>0)/length(maxcell)<0.5
        % turn matrix around
        maxcell=101-maxcell; dirdet=1;
    else dirdet=0;
    end
% make it one sequence - remove the restart
SIM= [t(maxtim) maxcell' maxval' ];
    [ SIMc, Cs ] = cont_slop( SIM(2:end,:), N,10000, ws, opt, mfreq); 
    if opt==1;
     figure; subplot(2,1,1); boxplot(Cs(:,2)); title('Distribution of Slope in Simulation');
             subplot(2,1,2); plot(Cs(:,1), Cs(:,2)/median(Cs(:,2)), 'b'); title('Slope and Frequency')
                 hold on; plot(Cs(:,1), Cs(:,5)/median(Cs(:,5)), 'k')
                 xlabel('Time'); ylabel('Normalized Slope and Frequency')
                 legend('Slope', 'Frequency', 'Location', 'Southeast')
           if length(unique(Cs(:,5)))>1;
                text(median(Cs(:,1)), min(Cs(:,5)/median(Cs(:,5))), ...
                    ['Freq from ', num2str(Cs(1,5)), ' to ', num2str(Cs(end,5))])
                subplot(2,1,1); boxplot(Cs(:,2),round(Cs(:,5)*100)./100); title('Distribution of Slope in Simulation');
           end
    end
 
end

