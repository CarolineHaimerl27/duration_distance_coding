function [ SIMc, C ] = cont_slop( SIM, N, Nseq, ws, opt_graph, mfreq )
% concatenate to one continuous slope from theta-modulated sequences
% assumes that cell indices is 1:N in the order in which they are presented
d1=1; d2=N; SIMc = SIM;
    for i=1:length(SIM(:,2))-1
        d1 = diff(SIMc(i:(i+1),2)); sd1 = -1*sign(d1); d2 = sd1.*(N-abs(d1));
        if abs(d1)>abs(d2);
            if sum(abs(d2)> abs((2:Nseq)*N-abs(d1)))>0; % abs((k+1)*N-abs(d1)); 
                indshub = find(min(abs((2:Nseq)*N-abs(d1)))==abs((2:Nseq)*N-abs(d1)));
                SIMc(i+1,2)= SIMc(i+1,2)+((1+indshub(1))*N*sd1);
            else
                if (SIMc(i+1,2)+(N*sd1))<0; 
                else SIMc(i+1,2)= SIMc(i+1,2)+(N*sd1);end
            end
            if abs(d2)==0;
                SIMc(i+1,2)= SIMc(i,2);
            end
            %if i<(length(SIM(:,2))-1); if (SIMc(i+1,2))>2*N & (SIM(i+1,2)-SIM(i+2,2))<0 ; k=k+1; end; end
            % k=k+1; at the point where the sequence switches for the x
            % time we need to asdjust!!!
        end
    end
if opt_graph==1; figure; plot(SIM(:,1), SIM(:,2), '-ob'); hold on; 
    plot(SIMc(:,1), SIMc(:,2), '-or'); end
% sliding-window regression
 step=1; if size(SIMc,1)<ws; ws=size(SIMc,1)-1; end
    for i=1:size(SIMc,1)-ws
        ci = robustfit(SIMc(i:(i+ws), 1), SIMc(i:(i+ws), 2));
        % save: time, coefficient, offset, starting cell
        C(i,:)=[median(SIMc(i:(i+ws), 1)) ci(2) ci(1) (SIM(i, 2)) median(mfreq(i:(i+ws))) SIMc(i, 1)];
    end
    
if opt_graph==1;   
    plot(SIMc(:,1), SIMc(:,1)*median(C(:,2))+min(C(:,3)), 'k');
    title('original slope, corrected slope, estimation'); xlabel('time(sec)'); ylabel('cell');
    legend('ori', 'corr', 'est'); end
end

% % make one continuous slope out of the sequence    
% d1=1; d2=N; SIMc = SIM;k=1;
%     for i=1:length(SIM(:,2))-1
%         d1 = diff(SIMc(i:(i+1),2)); sd1 = -1*sign(d1); d2 = sd1.*(k*N-abs(d1));
%         if abs(d1)>abs(d2);
%             if abs(d2)>abs((k+1)*N-abs(d1)); SIMc(i+1,2)= SIMc(i+1,2)+((k+1)*N*sd1);
%             else SIMc(i+1,2)= SIMc(i+1,2)+(k*N*sd1);end
%             if abs(d2)==0;
%                 SIMc(i+1,2)= SIMc(i,2);
%             end
%             %if i<(length(SIM(:,2))-1); if (SIMc(i+1,2))>2*N & (SIM(i+1,2)-SIM(i+2,2))<0 ; k=k+1; end; end
%             % k=k+1; at the point where the sequence switches for the x
%             % time we need to asdjust!!!
%         end
%     end
% figure; plot(SIM(:,1), SIM(:,2)); hold on; plot(SIMc(:,1), SIMc(:,2))
% % sliding-window regression
%     ws=10; step=1;
%     for i=1:size(SIMc,1)-10
%         ci = robustfit(SIMc(i:(i+10), 1), SIMc(i:(i+10), 2));
%         plot(SIMc(i:(i+10), 1), ci(1)+ci(2)*SIMc(i:(i+10), 1), 'k');
%         C(i,:)=[median(SIMc(i:(i+10), 1)) ci(2)];
%     end
%     title('original slope, corrected slope, estimation'); xlabel('time(sec)'); ylabel('cell');
%     legend('ori', 'corr', 'est')