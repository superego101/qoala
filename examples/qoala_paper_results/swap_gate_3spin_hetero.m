% Optimal control pulse optimisation for a swap gate over a scalar coupling 
% of three coupled 2-level systems, with each subsystem having its own 
% x-control and y-control. Controls are penalised with a spillout-norm
% square penalty.
%
%        140Hz         -160Hz
%   1H --------- 13C --------- 19F
%
% A comparison is made of wall-clock time and convergence, between of 
% GRAPE using QOALA and GRAPE using the exact Auxiliary Matrix method. Each
% of the different splitting methods are also run without adaptivity.
% Directional derivatives are calulated with a Taylor series version
% of the Auxiliary Matrix method.
% 
% Results are averaged from 28 different, random, initial pulse guess. The
% input "time-steps" below is an optional input, designed to be used when
% running these simulations over a number of different time slices.A good
% set of time slices to investigate are time_steps=[132 264 528 1056 2112].
%
% The reults of this simulation are using in FIG.4B and FIG.4D in the QOALA
% paper:                            http://dx.doi.org/......./.........
% 
% Auxiliary Matrix method:          http://dx.doi.org/10.1063/1.4928978
% Spillout-norm square penalty:     http://dx.doi.org/10.1063/1.4949534
%
% Author:       David L. Goodwin 
% Dev. period:  01/2020--04/2022
% Updated:      13/04/22
%
% Contact:
%   david.goodwin@chem.ox.ac.uk
%   david.goodwin@partner.kit.edu

function swap_gate_3spin_hetero(time_steps)

% simulate at these numbers of time slices
if nargin==1, time_steps_grid=time_steps;
elseif nargin==0, time_steps_grid=264;
end

% number of runs to average the optimisations
num_runs=28;

% data handling flags
plot_data=true;
save_data=true;

if plot_data
    if exist([cd filesep 'figures.zip'],'file'), unzip([cd filesep 'figures.zip']); delete([cd filesep 'figures.zip']); %#ok<*UNRCH> 
    elseif ~exist([cd filesep 'figures'],'dir'), mkdir([cd filesep 'figures'])
    end
end

if save_data
    if exist([cd filesep 'data.zip'],'file'), unzip([cd filesep 'data.zip']); delete([cd filesep 'data.zip']);
    elseif ~exist([cd filesep 'data'],'dir'),   mkdir([cd filesep 'data'])
    end
end

% bootstrap the inputs
omega=[0 0 0];          % offsets, Hertz
coupling=[140 -160];    % J-couplings, Hertz
pwr=1e3;                % Maximum pulse power, Hertz
T=26.4e-3;              % Total pulse duration, seconds
max_iter=250;           % maximum optimisation iterations

% add the paths
addpath(genpath('../../kernel/'))
addpath(genpath('../../utilities/'))

% grab the current file name
file_name = mfilename('fullpath'); [~,file_name]=fileparts(file_name); %#ok<ASGLU> 

% single spin-1/2 basis states (spherical tensors)
unit  = [+1  0  0  0].';
rho_z = [ 0  0 +1  0].';    rho_z=rho_z/2;

% single spin-1/2 Cartesian rotation operators
I  = [+1  0  0  0 ; 0 +1  0  0 ; 0  0 +1  0 ; 0  0  0 +1];
Jx = [ 0  0  0  0 ; 0  0 +1  0 ; 0 +1  0 +1 ; 0  0 +1  0];  Jx=Jx/sqrt(2);
Jy = [ 0  0  0  0 ; 0  0 -1i 0 ; 0 +1i 0 -1i; 0  0 +1i 0];  Jy=Jy/sqrt(2);
Jz = [ 0  0  0  0 ; 0 +1  0  0 ; 0  0  0  0 ; 0  0  0 -1];

% 2-spin rotation operators of composite system
LxH=sparse(kron(kron(Jx,I),I)); LyH=sparse(kron(kron(Jy,I),I)); LzH=sparse(kron(kron(Jz,I),I));
LxC=sparse(kron(kron(I,Jx),I)); LyC=sparse(kron(kron(I,Jy),I)); LzC=sparse(kron(kron(I,Jz),I));
LxF=sparse(kron(I,kron(I,Jx))); LyF=sparse(kron(I,kron(I,Jy))); LzF=sparse(kron(I,kron(I,Jz)));

% zz product operators
Lzz_12=kron(kron(Jz,unit*rho_z'+rho_z*unit'),I)+kron(kron(unit*rho_z'+rho_z*unit',Jz),I);
Lzz_13=kron(Jz,kron(I,unit*rho_z'+rho_z*unit'))+kron(unit*rho_z'+rho_z*unit',kron(I,Jz));
Lzz_23=kron(I,kron(Jz,unit*rho_z'+rho_z*unit'))+kron(I,kron(unit*rho_z'+rho_z*unit',Jz));

% contruct swap gate from Linden-1999 eq.14
U_swap= (expm(-1i*(+LzH+LzF)*3*pi/2)*...
        (expm(-1i*(-LxH)*pi/2)*...
        (expm(-1i*(+2*Lzz_13)*pi/2)*...
        (expm(-1i*(-LyH-LyF)*pi/2)*...
        (expm(-1i*(+2*Lzz_13)*pi/2)*...
        (expm(-1i*(+LxH+LxF)*pi/2)*...
        (expm(-1i*(+2*Lzz_13)*pi/2)*...
        (expm(-1i*(+LyH)*pi/2)))))))));

% clean to swap gate
U_swap=sparse(1e-9*round((1/1e-9)*(U_swap)));

% offsets in Zeeman Hamiltonian
H_offset=sparse(2*pi*omega(1)*LzH + 2*pi*omega(2)*LzC + 2*pi*omega(3)*LzF);

% J-coupling Hamiltonian
H_couple=sparse(2*pi*coupling(1)*Lzz_12+2*pi*coupling(2)*Lzz_23);

% Define control parameters
param.basis='sphten';                       % basis set
param.space='liouville';                    % state space
param.drift_sys.interaction=H_couple;       % Drift without offsets
param.operators={LxH,LyH,LxC,LyC,LxF,LyF};  % Control operators
param.prop_targ={U_swap};                   % Target propagator
param.pwr_levels=2*pi*pwr;                  % Pulse power (rad/s)
param.pulse_dur=T;                          % Pulse duration (s)
param.penalties={'SNS'};                    % Penalty method
param.prop_zeroed=1e-12;                    % tolerance to run expm
param.auxmat_method='taylor';               % taylor much faster than pade
param.nspins=length(omega);                 % number of spins in system

% define optimisation parameters
param.method='lbfgs';             % Optimisation method
param.n_grads=25;                 % lbfgs - gradient store number
param.tol_g=1e-12;                % Termination: norm(gradient,2)
param.tol_x=1e-12;                % Termination: norm(waveform change,1)
param.max_iter=max_iter;          % Termination: max iterations


% parameters for the auxiliary matrix method
sims_param{1}=param;
sims_param{1}.optimcon_fun=@grape_ugate_auxmat; % Optimal control function
sims_param{1}.drift_sys.singlespin=H_offset;    % Drift without coupling


% parameters for the qoala method
sims_param{2}=param;
sims_param{2}.prop_cache='carry';
sims_param{2}.optimcon_fun=@grape_ugate_qoala;      % Optimal control function
sims_param{2}.ctrl_axes={'x','y','x','y','x','y'};  % operation axis of controls
sims_param{2}.drift_sys.offset=2*pi*omega(:);       % Drift without coupling
sims_param{2}.nspins=length(omega);                 % number of spins
sims_param{2}.spin_control=[1 1 0 0 0 0;...         % ops 1,2 control 1st spin
                            0 0 1 1 0 0;...         % ops 3,4 control 2nd spin
                            0 0 0 0 1 1];           % ops 5,6 control 3rd spin
sims_param{2}.pauli_operators={LxH,LyH,LzH,[],[],[],[],[],[];
                               [],[],[],LxC,LyC,LzC,[],[],[];
                               [],[],[],[],[],[],LxF,LyF,LzF};

% parameters for the qoala method with fidelity check at each iteration
sims_param{2}.fidelity_chk=@waveform_fidelity;  % function to check fidelity
sims_param{2}.drift_sys.singlespin=H_offset;    % Drift without coupling

% splitting orders
sims_param{2}.drift_sys.split_order=2;
sims_param{3}=sims_param{2};
sims_param{3}.drift_sys.split_order=3;
sims_param{4}=sims_param{2};
sims_param{4}.drift_sys.split_order=4;
sims_param{5}=sims_param{2};
sims_param{5}.drift_sys.splitset=[2 3 4];

% save strings
nm={'auxmat','order2','order3','order4','qoala'}; %#ok<NASGU> 

% set up the figure
if plot_data
    for p=1:numel(sims_param)-1
        figure(p); clf;
        ln_colours=colororder({'#0072BD','#4DBEEE','#EDB120','#D95319','#A2142F'});
        if p==numel(sims_param)-1
            sgtitle('Three interacting 2-level systems with 6 controls, x and y on each system, adaptive')
        else
            sgtitle(['Three interacting 2-level systems with 6 controls, x and y on each system, split order = ' num2str(sims_param{p+1}.drift_sys.split_order)])
        end
        subplot(1,3,1)
        set(gca, 'YScale', 'log');
        xlim([0 max_iter]); ylim([1e-8 1]);
        xlabel('Iteraction'); ylabel('1-F'); title('Convergence');
        lgd=legend({},'location','northeast'); lgd.Title.String = 'QOALA'; grid minor; drawnow
        subplot(1,3,2)
        set(gca, 'YScale', 'log'); set(gca, 'XScale', 'log');
        xlim([0.02 10]); ylim([1e-8 1])
        xlabel('Wall-clock time (s)'); ylabel('1-F'); title('Computation')
        lgd=legend({},'location','southwest'); lgd.Title.String = 'Exact'; grid minor; drawnow
        subplot(1,3,3)
        set(gca, 'XScale', 'log');
        xlim([1e-5 1e-1]); ylim([1 25]); set(gca,'Xdir','reverse');
        xlabel('1-F'); ylabel('Speedup'); title('Comparative gain')
        legend({},'location','northeast'); grid minor; drawnow
    end
end

ctrl_struct=cell(1,numel(sims_param));

% grid of (in)fidelities to compare each method at these points
fidelity_check_space=logspace(-1,-9,91);

% initialise the results
inner_time=NaN*ones(numel(sims_param),param.max_iter+1,num_runs);
infidelity=NaN*ones(numel(sims_param),param.max_iter+1,num_runs);
fidelity_Q2=cell(1,numel(sims_param));      fidelity_Q2(:)={NaN*ones(param.max_iter+1,length(time_steps_grid))};
timer_Q2=cell(1,numel(sims_param));         timer_Q2(:)={ones(param.max_iter+1,length(time_steps_grid))};
timer_check_Q2=cell(1,numel(sims_param));   timer_check_Q2(:)={NaN*ones(length(fidelity_check_space),length(time_steps_grid))};

% run over all time steps set in the simulation grid
for ind=1:length(time_steps_grid)

    % parse the control parameters
    for p=1:numel(sims_param)
        sims_param{p}.pulse_nsteps=time_steps_grid(ind); % set the number of time steps
        ctrl_struct{p}=optimconset(sims_param{p}); % parse
    end

    % set num_runs different initial pulse guesses
    guess=cell(1,num_runs);
    for rndrun=1:num_runs
        rng(rndrun); % seed the random number generator
        guess{rndrun}=2*rand(numel(param.operators),time_steps_grid(ind)).'-1; % Initial guess
    end

    parfor rndrun=1:num_runs
        
        % initialise for the inner timers and fidelities
        timers=NaN*ones(numel(ctrl_struct),max_iter+1);
        fidelities=NaN*ones(numel(ctrl_struct),max_iter+1);

        % Run and time all the methods
        for p=1:numel(ctrl_struct)
            [~,data]=fmaxnewton(ctrl_struct{p},@optimfun_grape_xy,guess{rndrun});
            timers(p,:)=data.timer(:,4);                % pull the inner timer
            fidelities(p,:)=data.fx_chk_store(:,1);     % pull fidelity trajectory
        end
        inner_time(:,:,rndrun)=timers;          % store ther timers
        infidelity(:,:,rndrun)=1-fidelities;    % store the infidelities
    end

    % replace NaN of terminated optimisation with final fidelity achieved
    for p=1:numel(ctrl_struct), for k=1:num_runs, a=find(~isnan(squeeze(infidelity(p,:,k))),1,'last'); infidelity(p,a:end,k)=infidelity(p,a,k); end, end
    
    % ensure infidelity not less than eps
    infidelity(infidelity<eps)=eps;

    % calculate the medians
    fidelity_Q2{1}(:,ind)=prctile(squeeze(infidelity(1,:,:)),50,2);
    timer_Q2{1}(:,ind)=prctile(squeeze(inner_time(1,:,:)),50,2);

    for p=2:numel(ctrl_struct)
        fidelity_Q2{p}(:,ind)=prctile(squeeze(infidelity(p,:,:)),50,2);
        timer_Q2{p}(:,ind)=prctile(squeeze(inner_time(p,:,:)),50,2);

        for k=1:length(fidelity_check_space)
            a1=find(fidelity_Q2{p}(:,ind)<fidelity_check_space(k),1,'first');
            a2=find(fidelity_Q2{1}(:,ind)<fidelity_check_space(k),1,'first');
            if ~(isempty(a1)) && ~(isempty(a2))
                timer_check_Q2{p}(k,ind)=timer_Q2{p}(a1,ind);
                timer_check_Q2{1}(k,ind)=timer_Q2{1}(a2,ind);
            end
        end
    end

    if plot_data
        for p=2:numel(ctrl_struct)

            figure(p-1)
            subplot(1,3,1); hold on;
            semilogy(0:max_iter,fidelity_Q2{p}(:,ind),'Color',ln_colours(ind,:),'LineStyle','-', 'DisplayName', ['N=' num2str(time_steps_grid(ind))])
            semilogy(0:max_iter,fidelity_Q2{1}(:,ind),'Color',ln_colours(ind,:),'LineStyle','--', 'HandleVisibility','off')
            drawnow; hold off;

            subplot(1,3,2); hold on;
            loglog(timer_Q2{p}(:,ind),fidelity_Q2{p}(:,ind),'Color',ln_colours(ind,:),'LineStyle','-', 'HandleVisibility','off')
            loglog(timer_Q2{1}(:,ind),fidelity_Q2{1}(:,ind),'Color',ln_colours(ind,:),'LineStyle','--', 'DisplayName', ['N=' num2str(time_steps_grid(ind))])
            if (max(timer_Q2{1}(:,ind))>10) || (max(timer_Q2{p}(:,ind))>10), xlim([0.02 inf]); end
            drawnow; hold off;

            subplot(1,3,3); hold on;
            semilogx(fidelity_check_space,timer_check_Q2{1}(:,ind)./timer_check_Q2{p}(:,ind),'Color',ln_colours(ind,:),'DisplayName', ['N=' num2str(time_steps_grid(ind))],'LineWidth',2.0)
            if max(timer_check_Q2{1}(:,ind)./timer_check_Q2{p}(:,ind))>25, ylim([1 inf]); end
            drawnow; hold off;

        end
    end

    % re-evaluate for only the exact method
    for k=1:length(fidelity_check_space)
        a2=find(fidelity_Q2{1}(:,ind)<fidelity_check_space(k),1,'first');
        if ~(isempty(a2))
            timer_check_Q2{1}(k,ind)=timer_Q2{1}(a2,ind);
        end
    end

    if save_data

        for p=1:numel(ctrl_struct)

            % save escalade data to file
            conv_data=[(0:max_iter).' squeeze(inner_time(p,:,:)),...
                prctile(squeeze(inner_time(p,:,:)),25,2),...
                prctile(squeeze(inner_time(p,:,:)),50,2),...
                prctile(squeeze(inner_time(p,:,:)),75,2)];
            acc='   %12d'; headers={'iter'};
            for n=1:num_runs
                headers{n+1}=['rng(' int2str(n) ')']; acc=[acc '   %+12.5e']; %#ok<*AGROW>
            end
            headers{end+1}='Q1'; acc=[acc '   %+12.5e'];
            headers{end+1}='Q2'; acc=[acc '   %+12.5e'];
            headers{end+1}='Q3'; acc=[acc '   %+12.5e\n'];
            savetodat(['data' filesep file_name '_' num2str(time_steps_grid(ind),'%04.f') 'N_' nm{p} '_time.dat'],headers,(conv_data),acc);

            % save escalade check data to file
            conv_data=[(0:max_iter).' squeeze(infidelity(p,:,:)),...
                prctile(squeeze(infidelity(p,:,:)),25,2),...
                prctile(squeeze(infidelity(p,:,:)),50,2),...
                prctile(squeeze(infidelity(p,:,:)),75,2)];
            acc='   %12d'; headers={'iter'};
            for n=1:num_runs
                headers{n+1}=['rng(' int2str(n) ')']; acc=[acc '   %+12.5e'];
            end
            headers{end+1}='Q1'; acc=[acc '   %+12.5e'];
            headers{end+1}='Q2'; acc=[acc '   %+12.5e'];
            headers{end+1}='Q3'; acc=[acc '   %+12.5e\n'];
            savetodat(['data' filesep file_name '_' num2str(time_steps_grid(ind),'%04.f') 'N_' nm{p} '_infidelity.dat'],headers,(conv_data),acc);
        end
    end

    % inform the user
    disp(['Finished optimisations for ' num2str(time_steps_grid(ind)) ' time steps'])

end

if plot_data
    for p=2:numel(ctrl_struct)
        figure(p-1)
        subplot(1,3,3); hold on;
        semilogx(fidelity_check_space,(min(timer_check_Q2{1},[],2)./min(timer_check_Q2{p},[],2)),'Color',[0 0 0],'LineStyle',':','DisplayName', 'all','LineWidth',2.0)
        if max((min(timer_check_Q2{1},[],2)./min(timer_check_Q2{p},[],2)))>25, ylim([1 inf]); end
        drawnow; hold off;
        savefig(['figures' filesep file_name '_' num2str(time_steps_grid(ind),'%04.f') 'N_' nm{p}]);
    end
    zip([cd filesep 'figures.zip'],[cd filesep 'figures']); rmdir([cd filesep 'figures'],'s');
end
if save_data
    zip([cd filesep 'data.zip'],[cd filesep 'data']); rmdir([cd filesep 'data'],'s');
end

disp('FIN')

end


% cost function
function [trajdat,fdty,grad]=optimfun_grape_xy(wf,param)

% Count the outputs
n_outputs=nargout;

% check an ensemble is not present
if numel(param.drift_sys)==1, drift_terms=param.drift_sys{1};
else, error('this function is not coded for an ensemble')
end

% hack to be fixed in optimconset
if contains(func2str(param.optimcon_fun),{'escalade','qoala'})
    ctrls=param.pauli_operators;
    optimfun=param.optimcon_fun;
else
    ctrls=param.operators;
    optimfun=param.optimcon_fun;
end

% extract the power levels and fidelity type
nsteps=size(wf,1); kctrls=size(wf,2);

% Get initial and target state
init_state=param.initials{1};
targ_state=param.targets{1};

% Get the power_level
power_lvls=param.pwr_levels(1,:);

% Count penalty terms
npenterms=numel(param.penalties);

% Decide how much needs computing
if n_outputs==2

    % Preallocate output
    fdty=zeros(1,npenterms+1);

    % Calculate the objective and its gradient
    [trajdat,fdty(1)]=optimfun(param,drift_terms,ctrls,power_lvls,wf,init_state,targ_state);

    % Apply all penalties
    for n=1:npenterms

        % Calculate the penalty
        fdty(n+1)=penalty_fun(param,wf,param.penalties{n},param.l_bound,param.u_bound,param.p_weights(n));

    end

elseif n_outputs==3

    % Preallocate output
    fdty=zeros(1,npenterms+1);
    grad=zeros(nsteps,kctrls,npenterms+1);

    % Calculate the objective and its gradient
    [trajdat,fdty(1),grad(:,:,1)]=optimfun(param,drift_terms,ctrls,power_lvls,wf,init_state,targ_state);

    % Apply all penalties
    for n=1:npenterms

        % Calculate the penalty
        [fdty(n+1),grad(:,:,n+1)]=penalty_fun(param,wf,param.penalties{n},param.l_bound,param.u_bound,param.p_weights(n));

    end

end

% small hack, but may be useful for an eventual ensemble
if isfield(trajdat,'drift_sys'), trajdat.drift_sys={trajdat.drift_sys}; end

end

function savetodat(filename,headers,datacols,acc) %#ok<DEFNU> 

% Write to file
fileID = fopen(filename,'w');
for n=1:numel(headers)
    fprintf(fileID, pad(headers{n},15,'left'));
end
fprintf(fileID, '\n');
for n=1:size(datacols,1)
    fprintf(fileID, acc, datacols(n,:));
end

fclose(fileID);

end
