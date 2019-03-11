classdef FRModel < handle
    
    properties (SetAccess = private)
        N
        inpc
        tend
        dt
        Ufix
        tauR
        tauF
        Mscale
        thetaon
        freq
        inpth
        
        dirdet
        M
        Udyn
        Xdyn
        t
        slope
        thslope
        medthslope
        thshift
        C
        maxfir 
        
    end
    
    methods
                
        function this = FRModel(inpc, tend, dt, Mscale, thetaon, freq, inpth, tbreak, modelspec)
            
            %%%%%%%%%% INPUTS CHECKS %%%%%%%%%%%%
            mfname = mfilename;

            if nargin < 1 || isempty(inpc)
                warning('%s: the user must input a constant input', mfname);
            end
            if nargin < 2 || isempty(tend)
                warning('%s: the user must input the length of the simulation in sec', mfname);
            end
            if nargin < 3 || isempty(dt)
                warning('%s: no dt specified, use default 1e-4', mfname);
                dt=1e-4;
            end
            if nargin < 4 || isempty(Mscale)
                warning('%s: no scale for starting firing rate, use default 0.1', mfname);
                Mscale=0.1;
            end
            if nargin < 5 || isempty(thetaon)
                warning('%s: theta input set to default, off', mfname);
                thetaon=0; freq=0; inpth=0;
            end
            if nargin < 6 || isempty(freq)
                warning('%s: no theta frequency specified, theta input put off', mfname);
                thetaon=0; freq=0; inpth=0;
            end
            if nargin < 7 || isempty(inpth)
                warning('%s: no theta amplitude specified, theta input put off', mfname);
                thetaon=0; freq=0; inpth=0;
            end
            if nargin < 8 || isempty(tbreak)
                warning('%s: no break time specified, set to 0', mfname);
                tbreak=0;
            end
            if thetaon ==0; freq=0; inpth=0; end
            global N Ufix tauR tauF
            X0 = ones(N,1); U0 = ones(N,1)*Ufix; M0=[linspace(0, 1, N)*Mscale]'; %rand(N,1)
            % copy input parameters to object properties
            this.N      = N;
            this.inpc   =inpc;
            this.tend   = tend;
            this.dt     = dt;
            this.Ufix   = Ufix;
            this.tauR   = tauR;
            this.tauF   = tauF;
            this.Mscale = Mscale;
            this.thetaon = thetaon;
            this.freq   = freq;
            this.inpth  = inpth;
            % define parameters for oscillatory input
            I_theta=inpth;
            f_theta = freq;
            par = [0, inpc, freq, 0, inpth];
            % compute firing rate dynamics
            opts = odeset('Events',@Events);
            model = str2func(modelspec);
            [t1,x]=ode45(model,0:dt:tend,[par'; X0; U0; M0], opts); % Runge Kutta integration
            M = x(:,length(par)+1+2*N:end);
            this.Udyn = x(:,length(par)+1+N:end-N);
            this.Xdyn = x(:,length(par)+1:end-2*N);
            % determine break behavior
            if tbreak > 0
                % simulation of break
                X       = this.Xdyn(end,:);%x(:,length(par)+1:length(par)+N);
                Uendb = this.Udyn(end,:);
                par([2,5])  = 0;
                [t2,x]  = ode45(model,tend:dt:(tend+tbreak),[par'; X'; Uendb'; M(end, :)'], opts); % Runge Kutta integration
                Mbreak  = x(:,length(par)+1+2*N:end);
                Ubreak = x(:,length(par)+1+N:end-N);
                Xbreak  = x(:,length(par)+1:length(par)+N);
                this.Udyn = [this.Udyn; Ubreak];
                this.Xdyn = [this.Xdyn; Xbreak];
                % simulation after break
                par([2,5])  = [inpc, I_theta];
                [t3,x]  = ode45(model,(tend+tbreak):dt:(2*tend+tbreak),[par'; Xbreak(end, :)'; Ubreak(end, :)'; Mbreak(end, :)'], opts); % Runge Kutta integration
                Mabreak = x(:,length(par)+1+2*N:end);
                %                         % determine slope after break
                %                         [ ~, Cs, ~, this.dirdet(ia, ifr,2) ] = SLOPE( Mabreak, t3, f_theta, I_theta, N, size(Mabreak,1)/2,0);
                %                         this.C(ia, ifr, 2)=Cs(2);
                this.M = [M; Mbreak; Mabreak]; this.t = [t1; t2; t3];
                this.Udyn = [this.Udyn; x(:,length(par)+1+N:end-N)];
                this.Xdyn = [this.Xdyn; x(:,length(par)+1:length(par)+N)];
            else
                this.M = M; this.t = t1;
            end
            
            % Average Run Sequence Slope
            % if there was a break - compute the slope before and after the break
            if tbreak == 0
                Mtmp = M;%(:,1:100);
                ttemp = 0:dt:range(t1);
                % N1 = 100;
            else
                Mtmp = Mabreak;
                ttemp = 0:dt:range(t3);
                % N1 = 100;
            end
            this.slope(1)                  = RunSeqSlope(Mtmp, N, ttemp', thetaon, f_theta);
            [ ~, Cs, ~, this.dirdet ]   = SLOPE( Mtmp, ttemp', f_theta, I_theta, N, size(Mtmp,1)/2,0);
            this.C(1)=Cs(2);
            % Theta slope
            [this.thslope, this.thshift]= ThSlope(ttemp, f_theta, N, Mtmp, dt);
            this.medthslope(1)             = median(this.thslope, 'omitnan');
            % compute median max firing rate
            this.maxfir(1)                 = median(max(Mtmp));
            if tbreak >0
                    Mtmp2 = M;
                    ttemp = 0:dt:range(t1);
                    this.slope(2)                  = RunSeqSlope(Mtmp2, N, ttemp', thetaon, f_theta);
                    [ ~, Cs]   = SLOPE( Mtmp2, ttemp', f_theta, I_theta, N, size(Mtmp2,1)/2,0);
                    this.C(2) = Cs(2);
                    % Theta slope
                    x= ThSlope(ttemp, f_theta, N, Mtmp2, dt);
                    this.medthslope(2)             = median(x, 'omitnan');
                    % compute median max firing rate
                    this.maxfir(2)                 = median(max(Mtmp2));
            end
                
            
        end
    end
    
    
    
end
