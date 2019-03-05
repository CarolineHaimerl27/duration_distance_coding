clear all; close all; clc

transient = 1;

% Grid Size

MM = 100;

% Number of neurons
N = 50;

% Parameters

J1 = 70; % strength of map-specific interaction (excitatory strength)
J0 = 50; % uniform feedback inhibition (inhibitory strength)
Ufix=0.15; 
tauF = 0.3; %facilitation recovery time constant
tauR = 0.3; % Depression time constant

% Frequency of the Theta forcing

fTheta = 9; % In Herz

% Membrane time constant
tau = 0.01;

% Initialize the matrix where the data will be stored
DataMatrix = zeros(MM);

% Minimum Time required to find the peaks is the time taken by 40 cycles

NumTheta_Cycles = 40;
time_peaks = NumTheta_Cycles/fTheta;
t_minimum = transient + time_peaks;

% Grid in delta / I_theta parameter space
deltaVec = linspace(0.01,0.05,MM);
InputVec = linspace(0.1,20,MM);

for jj = 1:MM
    tic
    delta = deltaVec(jj);
    parfor ll = 1:MM        
        
        % Intial conditions
        U0 = Ufix*ones(N,1); M0=rand(N,1); D0 = rand(N,1);
        Ipar = InputVec(ll);
        
        parameters = [N J1 J0 tauF tauR tau Ufix Ipar fTheta delta];
        X0 = [M0;U0;D0];
        [time_sim,x_sim] = Basic_Model(parameters,X0,[0 t_minimum]);
        
        if(~isempty(x_sim))
            
            % Eliminate the first part of Transient
            
            oo = find(time_sim>transient);
            time_sim = time_sim(oo);
            M = x_sim(oo,1:N);
            
            % Select the positive parts of the theta sequence in which, by
            % definition I will have the maximum activity.
            
            Input = sin(2*pi*fTheta*time_sim);
            oo = Input<0;
            
            % Identify the time indeces of changes from positive to negative
            % and viceversa:
            
            aux = diff(oo);
            negativeToPositive=(find(aux==-1));
            negativeToPositive=negativeToPositive+1;
            
            positiveToNegative=(find(aux==1));
            positiveToNegative=positiveToNegative+1;
            
            % Make sure that we have as many positive->negative changes as
            % negative -> positive.
            
            if(positiveToNegative(1)<negativeToPositive(1))
                positiveToNegative(1)=[];
            end
            
            t_max = zeros(1,NumTheta_Cycles-1);
            neu_max = zeros(1,NumTheta_Cycles-1);
            
            for kk=1:NumTheta_Cycles-1
                
                % For each positive theta cycle I choose the neuron with maximum
                % value of the together with the time it happened
                
                time_chunk = time_sim(negativeToPositive(kk):positiveToNegative(kk));
                Data_chunk = M(negativeToPositive(kk):positiveToNegative(kk),:);
                [~,IndMax] = max(Data_chunk(:));
                [I_time, I_Neuron] = ind2sub(size(Data_chunk),IndMax);
                t_max(kk) = time_chunk(I_time); neu_max(kk) = I_Neuron;
                
            end
            
            % Make sure that the periodic boundary conditions are met:
            
            zz = diff(neu_max);
            hh = find(zz<0);
            while(~isempty(hh))
                neu_max(hh(1)+1:end) = neu_max(hh(1)+1:end)+N;
                zz = diff(neu_max);
                hh = find(zz<0);
            end
            
            % Calculate the propagation velocity of the bump for the value
            % of delta and input evaluated in the loop
            DataMatrix(jj,ll) = mean(2*pi/N*diff(neu_max)./diff(t_max));
            
        else
            
            % If something went wrong with simulation I set the value of
            % propagation to 0.
            DataMatrix(jj,ll) = 0;
        end
                
        
    end
    toc
end


%%

figure
hold all
for i = 1:MM
    plot(InputVec,DataMatrix(i,:))
end
    
figure
surf(InputVec,deltaVec,DataMatrix')

% Save the Data MAtrix in PhaseDiagram
save('PhaseDiagram.mat','DataMatrix','InputVec','deltaVec');