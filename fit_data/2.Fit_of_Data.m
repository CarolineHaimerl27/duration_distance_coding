clear all; close all; clc

% Load the data that we run previously on main_model
load('PhaseDiagram.mat');

% Plot the 2D diagram
figure
surf(InputVec,deltaVec,MatrixOfStuff')
shading flat
view(2)
ylim([0.005 0.1])
xlim([0.1 3])
set(gca,'CLim',[0 4])

% Initialize the alpha and beta parameters I_theta = alpha + beta*speed

alph = linspace(0.01,10,200);
bet = linspace(0.001,1,200);

% This is the vector with the recorded data
VectorData =   [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% For each recording session this is the result of the coding test (1-5). 
Nature =       [1 2 3 3 2 5 4 0 1  0  0  3  3  2  4  2  2  3  3  4  5  4  5  4  5  4  4  5  2 2  3  3  2  3];
% Number of Cells recorded on each session
Ncell =        [171 74 71 65 66 36 31 38 30 32 23 22 32 22 33 40 23 52 37 63 65 56 54 55 57 32 29 70 109 72 79 114 69 61];

% Initialize the matrix were the fitted (alpha, beta and delta) will be stored
Parameters = zeros(length(VectorData),3); %

figure
counter = 1;

% This is the minimization routine that fits the parameter for each
% recording session

for i = VectorData
    
    load(['Data/' num2str(i) '.mat']);
    
    [SpeedSort,idx] = sort(Speed);
    VelSort = Velocity(idx)*2*pi/Ncell(i);
    
    ErrorMin = 10000*ones(length(deltaVec),1);
    alMin = 1000*ones(length(deltaVec),1);
    beMin = 1000*ones(length(deltaVec),1);
    
    tic
    for ll = 1:1
        
        meanValueIdx = round(size(MatrixOfStuff,1)/10);
        
        EE = 100000;
        aa = Inf;
        bb = Inf;
        VelSim = MatrixOfStuff(meanValueIdx,:);
        
        for j = 1:length(alph)
            for k = 1:length(bet)
                
                I = alph(j) + bet(k)*SpeedSort;
                if(max(I)<15 && min(I)>0.1)
                    iid = nearestpoint(I,InputVec);
                    Err = sum(abs((VelSort' - VelSim(iid))));
                    
                    
                    
                    if(Err<EE)
                        EE = Err;
                        aa = alph(j);
                        bb = bet(k);
                    end
                    
                else
                    break
                end
                
            end
        end
        
        ErrorMin(ll) = EE;
        alMin(ll) = aa;
        beMin(ll) = bb;
        
    end
    toc
    
    idx2 = find(ErrorMin==min(ErrorMin));
    
    Parameters(counter,1) = deltaVec(idx2);
    Parameters(counter,2) = alMin(idx2);
    Parameters(counter,3) = beMin(idx2);
    
    Idef = alMin(idx2) + beMin(idx2)*SpeedSort;
    
    subplot(7,6,counter)
    hold all
    plot(InputVec,MatrixOfStuff(meanValueIdx,:))
    plot(Idef,VelSort,'or')
    %xlim([0.1 3])
    
    counter = counter +1;
    medianSpeed(i) = median(SpeedSort);
    
end


% Save the data in Parameters.mat

save('Parameters.mat','Parameters','VectorData')

