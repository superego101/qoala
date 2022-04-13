% Gradient Ascent Pulse Engineering (GRAPE) objective function and 
% gradient. Adaptively set the Trotter number and splitting order, for a
% solver based split ESCALADE method.
%
% Inputs:
%
%       param       - general parameters
%
%       drift_sys   - information of the drift
%
%       ctrls       - control opterators
%
%       amps        - amplitudes of each control channel (rad/s)
%
%       wf          - piecewise constant waveform
%
%       init        - initial propagator
%
%       targ        - target effective propagator
%
% Outputs:
%
%       data        - structure containing data passed back to optimiser
%
%       fdty        - fidelity
%
%       grad        - gradient of fidelity wrt waveform controls
%
% Author:       David L. Goodwin 
% Dev. period:  01/2020--04/2022
% Updated:      13/04/22
%
% Contact:
%   david.goodwin@chem.ox.ac.uk
%   david.goodwin@partner.kit.edu

function  [data,fdty,grad]=grape_ugate_qoala(param,drift_sys,ctrls,amps,wf,init,targ)

% number of outputs
n_outputs=nargout;

% check if single spin simulation
if ~isfield(drift_sys,'interaction')
    if n_outputs==2,        [data,fdty]=escalade(param,drift_sys,ctrls,amps,wf,init,targ);
    elseif n_outputs==3,    [data,fdty,grad]=escalade(param,drift_sys,ctrls,amps,wf,init,targ);
    end; return;
end

% start the inner timer
esc_time=tic;

% Decide how much needs computing
if n_outputs==2

    % shorthands
    p_order=drift_sys.split_order;  q_trott=drift_sys.trotter_number;

    if strcmp(param.prop_cache,'carry') && isfield(drift_sys,'adaptset')
        
        % pull the adaptive bounds
        adapt_set=drift_sys.adaptset;

        % find where we are on the adaptivity variables
        k=find(ismember(adapt_set.',[q_trott; p_order].','rows'));

        drift_sys_in=drift_sys;
        drift_sys_in.spn_mults=drift_sys_in.spn_mults{k};
        drift_sys_in.ind_spn=drift_sys_in.ind_spn{k};
        drift_sys_in.ind_int=drift_sys_in.ind_int{k};
        drift_sys_in.inter_prop=drift_sys_in.inter_prop(k,:);
        drift_sys_in.inter_prop=drift_sys_in.inter_prop(~cellfun('isempty',drift_sys_in.inter_prop));
    else
        drift_sys_in=drift_sys;
    end

    fdty=qoala_fidelity(param,drift_sys_in,ctrls,amps,wf,init,targ);

elseif n_outputs==3 && ~drift_sys.adapt

    % shorthands
    p_order=drift_sys.split_order;  q_trott=drift_sys.trotter_number;

    if strcmp(param.prop_cache,'carry') && isfield(drift_sys,'adaptset')
        
        % pull the adaptive bounds
        adapt_set=drift_sys.adaptset;

        % find where we are on the adaptivity variables
        k=find(ismember(adapt_set.',[q_trott; p_order].','rows'));

        drift_sys_in=drift_sys;
        drift_sys_in.spn_mults=drift_sys_in.spn_mults{k};
        drift_sys_in.ind_spn=drift_sys_in.ind_spn{k};
        drift_sys_in.ind_int=drift_sys_in.ind_int{k};
        drift_sys_in.inter_prop=drift_sys_in.inter_prop(k,:);
        drift_sys_in.inter_prop=drift_sys_in.inter_prop(~cellfun('isempty',drift_sys_in.inter_prop));
    else
        drift_sys_in=drift_sys;
    end

    [data,fdty,grad]=qoala_gradient(param,drift_sys_in,ctrls,amps,wf,init,targ);

elseif n_outputs==3

    % if we are set as adaptive, but the number of iteractions since the
    % last adaptive check is less than adapt_minit, run escalade without
    % changing the splitting parameters.
    if drift_sys.adapt_counter < drift_sys.adapt_minit

        % shorthands
        p_order=drift_sys.split_order;  q_trott=drift_sys.trotter_number;

        % pull the adaptive bounds
        adapt_set=drift_sys.adaptset;

        if strcmp(param.prop_cache,'carry')

            % find where we are on the adaptivity variables
            k=find(ismember(adapt_set.',[q_trott; p_order].','rows'));

            drift_sys_in=drift_sys;
            drift_sys_in.spn_mults=drift_sys_in.spn_mults{k};
            drift_sys_in.ind_spn=drift_sys_in.ind_spn{k};
            drift_sys_in.ind_int=drift_sys_in.ind_int{k};
            drift_sys_in.inter_prop=drift_sys_in.inter_prop(k,:);
            drift_sys_in.inter_prop=drift_sys_in.inter_prop(~cellfun('isempty',drift_sys_in.inter_prop));
        else
            drift_sys_in=drift_sys;
        end

        % turn off the adaptivity flag
        drift_sys_in.adapt=false;

        % calculate fidelity and gradient
        [data,fdty,grad]=qoala_gradient(param,drift_sys_in,ctrls,amps,wf,init,targ);

        data.drift_sys=drift_sys;

        % increment the non-adaptive iteraction counter
        data.drift_sys.adapt_counter=drift_sys.adapt_counter+1;

        % turn on the adaptivity flag
        data.drift_sys.adapt=true;

    else
        % pull the size of the Hamiltonian
        dim=param.dim;
        
        % function counter
        f_adapt_count=0;

        % shorthands
        p_order=drift_sys.split_order;  q_trott=drift_sys.trotter_number;
        adapt_method=drift_sys.adapt_method;

        % pull the adaptive bounds
        adapt_set=drift_sys.adaptset;

        if ismember(adapt_method,{'exact'})

            % run an exact fidelity check to get the final state
            [data_exact,~]=waveform_fidelity(param,drift_sys,ctrls,amps,wf,init,targ);

            % pull the final state/propagator
            rho_0=data_exact.trajectory.final; %rho_0=rho_0/norm(full(rho_0),2);

        end

        if strcmp(param.prop_cache,'carry')

            % find where we are on the adaptivity variables
            k=find(ismember(adapt_set.',[q_trott; p_order].','rows'));

            drift_sys_in=drift_sys;
            drift_sys_in.spn_mults=drift_sys_in.spn_mults{k};
            drift_sys_in.ind_spn=drift_sys_in.ind_spn{k};
            drift_sys_in.ind_int=drift_sys_in.ind_int{k};
            drift_sys_in.inter_prop=drift_sys_in.inter_prop(k,:);
            drift_sys_in.inter_prop=drift_sys_in.inter_prop(~cellfun('isempty',drift_sys_in.inter_prop));
        else
            drift_sys_in=drift_sys;
        end

        % initial run to calculate the objective (don't add to counter)
        [fidelity_est,fwd_traj,fwd_split,P_n]=qoala_fidelity(param,drift_sys_in,ctrls,amps,wf,init,targ);

        % store the old trajectories and propagators
        P_str=P_n;  fwd_str=fwd_traj;   splt_str=fwd_split;

        % pull the final state and normalise
        rho_1=fwd_traj{end};
        rho_est=rho_1;%/norm(full(rho_1),2);

        % initialise the error
        error_f=1;

        % calculate the error for an exact evaluation
        if ismember(adapt_method,{'exact'}), error_f=error_f*abs(1-real(trace(rho_est'*rho_0)/dim)); end

        % set the threshold for adaptivity
        tol_f=drift_sys.adapt_tols(1);

        % decide to scale the adaptivity tolerance
        if ismember({'gradient'},drift_sys.adapt_scale) && isfield(param,'grad_norm')

            % scale with the norm or the last gradient
            tol_f=min(param.grad_norm,1)*tol_f;

        end
        if ismember({'infidelity'},drift_sys.adapt_scale)

            % scale with the fidelity estimate
            tol_f=max(eps,abs(1-max(0,min(fidelity_est,1))))*tol_f;
        end

        % find where we are on the adaptivity variables
        k=find(ismember(adapt_set.',[q_trott; p_order].','rows'));

        % loop until error is less than the tolerance
        while (error_f > tol_f)

            % store the old trajectories and propagators
            P_str=P_n;  fwd_str=fwd_traj;   splt_str=fwd_split;

            % increment the index and break if beyond bounds
            k=k+1; if k>size(adapt_set,2), break; end

            % use the lats final state
            if ismember(adapt_method,{'gain'}), rho_0=rho_1; end

            % next Trotter number and split order
            q_trott=adapt_set(1,k);     p_order=adapt_set(2,k);

            % find the constants and indices
            if ~strcmp(param.prop_cache,'carry')
                param.prop_method='taylor';
                drift_sys_in=drift_sys;
                [spn_c,ind_spn,ind_int,prop]=operator_splittings(param,p_order,q_trott,param.pulse_dt(1),drift_sys_in.interaction);

                % modify the split parameters
                drift_sys_in.split_order=p_order;    drift_sys_in.spn_mults=spn_c;
                drift_sys_in.ind_spn=ind_spn;            drift_sys_in.ind_int=ind_int;
                drift_sys_in.trotter_number=q_trott;
                drift_sys_in.inter_prop=prop;
            else

                drift_sys_in=drift_sys;
                drift_sys_in.spn_mults=drift_sys_in.spn_mults{k};
                drift_sys_in.ind_spn=drift_sys_in.ind_spn{k};
                drift_sys_in.ind_int=drift_sys_in.ind_int{k};
                drift_sys_in.inter_prop=drift_sys_in.inter_prop(k,:);
                drift_sys_in.inter_prop=drift_sys_in.inter_prop(~cellfun('isempty',drift_sys_in.inter_prop));
            end

            % calculate the final state/propagator
            [fidelity_est,fwd_traj,fwd_split,P_n]=qoala_fidelity(param,drift_sys_in,ctrls,amps,wf,init,targ);

            % pull the final state
            rho_1=fwd_traj{end};
            rho_est=rho_1;%/norm(full(rho_1),2);

            % increment the counter
            f_adapt_count=f_adapt_count+1;

            % calculate the overlap with previous final state/propagator
            error_f=abs(1-real(trace(rho_est'*rho_0)/dim));

            % reset the threshold for adaptivity
            tol_f=drift_sys.adapt_tols(1);

            % decide to scale the adaptivity tolerance
            if ismember({'gradient'},drift_sys.adapt_scale) && isfield(param,'grad_norm')

                % scale with the norm or the last gradient
                tol_f=min(param.grad_norm,1)*tol_f;

            end

            if ismember({'infidelity'},drift_sys.adapt_scale)

                % scale with the fidelity estimate
                tol_f=max(eps,abs(1-max(0,min(fidelity_est,1))))*tol_f;
            end

            % counter was set to inf before first iteration -- allow many
            % jumps at the first iteration to avoid an infeasible start.
            if ~isinf(drift_sys.adapt_counter)

                % little interval bug -- different for 'gain' and 'exact'
                if ismember(adapt_method,{'gain'}),         if f_adapt_count>drift_sys.adapt_interval, break; end
                elseif ismember(adapt_method,{'exact'}),    if f_adapt_count>drift_sys.adapt_interval-1, break; end
                end
            end
        end

        if ismember(adapt_method,{'exact'})

        else
            % previous Trotter number and split order
            if k>1
                q_trott=adapt_set(1,k-1);	p_order=adapt_set(2,k-1);
                if ~strcmp(param.prop_cache,'carry')
                    param.prop_method='taylor';
                    [spn_c,ind_spn,ind_int,prop]=operator_splittings(param,p_order,q_trott,param.pulse_dt(1),drift_sys.interaction);

                    % modify the split parameters
                    drift_sys.split_order=p_order;    drift_sys.spn_mults=spn_c;
                    drift_sys.ind_spn=ind_spn;            drift_sys.ind_int=ind_int;
                    drift_sys.trotter_number=q_trott;
                    drift_sys.inter_prop=prop;

                    % reset the split counter
                    drift_sys.adapt_counter=1;

                    % turn off adaptive section
                    drift_sys.adapt=false;
                    drift_sys_in=drift_sys;
                else

                    % modify the split parameters
                    drift_sys.split_order=p_order;
                    drift_sys.trotter_number=q_trott;

                    % reset the split counter
                    drift_sys.adapt_counter=1;

                    % turn off adaptive section
                    drift_sys.adapt=false;
                    drift_sys_in=drift_sys;
                    drift_sys_in.spn_mults=drift_sys_in.spn_mults{k-1};
                    drift_sys_in.ind_spn=drift_sys_in.ind_spn{k-1};
                    drift_sys_in.ind_int=drift_sys_in.ind_int{k-1};
                    drift_sys_in.inter_prop=drift_sys_in.inter_prop(k-1,:);
                    drift_sys_in.inter_prop=drift_sys_in.inter_prop(~cellfun('isempty',drift_sys_in.inter_prop));
                end
            end
        end

        % fidelity and gradient (store in param.adaptive)
        [data,fdty,grad]=qoala_gradient(param,drift_sys_in,ctrls,amps,wf,init,targ,fwd_str,splt_str,P_str);

        % turn on adaptive section
        drift_sys.adapt=true;

        % add the modified param to the data structure
        data.drift_sys=drift_sys;

        % store the counter
        data.f_adapt_count=f_adapt_count;

    end

    % store the grad norm
    data.grad_norm=norm(grad,2);

end

% stop the timer
data.timer=toc(esc_time);

end

% create a propagator matrix given indices and non-zero elements
function P_out=ind2propagator(prop_ind,P_in,n,nspins,unit,dim,method,dense_matrix)

% pull the required row for this propagator
P_n=P_in(n,:);

% all elements, including zeros, stored -- form the matrix
if numel(P_n)==dim*dim

    % form matrix from vector
    P_out=reshape(P_n,dim,dim);

    % decide if full or sparse
    if nnz(P_out)>dim*dim*dense_matrix, P_out=full(P_out); end
    return;
end

if strcmp(method,'prod')
    if numel(P_n)==size(prop_ind,1)-numel(unit)
        P_out=sparse(prop_ind(:,1),prop_ind(:,2),[unit; P_n(:)]);
    elseif size(P_n,2)==size(prop_ind,1)-numel(unit)
        P_out=1;
        for k=1:nspins
            P_1=P_in(n(k),:);
            P_1=sparse(prop_ind(:,1),prop_ind(:,2),[unit; P_1(:)]);
            P_out=(kron(P_out,P_1));
        end
    else
        error('incorrect number of propagator indices')
    end
elseif strcmp(method,'sum')
    prop_ind=cell2mat(prop_ind.');
    if numel(P_n)==size(prop_ind,1)-numel(unit)
        P_out=sparse(prop_ind(:,1),prop_ind(:,2),[unit; P_n(:)]);
    elseif numel(P_n)*4^(nspins-1)==size(prop_ind,1)
        P_1=repmat(P_in(:,n), 4^(nspins-1),1);
        P_out=sparse(prop_ind(:,1),prop_ind(:,2),P_1(:));
    elseif size(P_n,1)==size(prop_ind,1)
        P_out=spalloc(dim,dim,size(P_in,1)*(dim^(nspins-1)));
        for k=1:nspins
            P_2=1;
            for j=1:nspins
                if j==k
                    P_1=P_in(:,n(k));
                    P_1=sparse(prop_ind(:,1),prop_ind(:,2),[unit; P_1(:)]);
                else
                    P_1=speye(dim^(1/nspins));
                end
                P_2=kron(P_2,P_1);
            end
            P_out=P_out+P_2;
        end
    else
        error('incorrect number of propagator indices')
    end
end

% decided if full or sparse
if nnz(P_out)>dim*dim*dense_matrix, P_out=full(P_out); end

end

% create a propagator matrix given indices and non-zero elements
function P_out=ind2propagatorPROP(prop_ind,P_in,n,nspins,unit,dim,method,dense_matrix)

% pull the required row for this propagator
P_n=P_in(:,n);

% all elements, including zeros, stored -- form the matrix
if numel(P_n)==dim*dim

    % form matrix from vector
    P_out=reshape(P_n,dim,dim);

    % decide if full or sparse
    if nnz(P_out)>dim*dim*dense_matrix, P_out=full(P_out); end
    return;
end

if strcmp(method,'prod')
    if numel(P_n)==size(prop_ind,1)-numel(unit)
        P_out=sparse(prop_ind(:,1),prop_ind(:,2),[unit; P_n(:)]);
    elseif size(P_n,1)==size(prop_ind,1)
        P_out=1;
        for k=1:nspins
            P_1=P_in(:,n(k));
            P_1=sparse(prop_ind(:,1),prop_ind(:,2),P_1);
            P_out=(kron(P_out,P_1));
        end
    else
        error('incorrect number of propagator indices')
    end
elseif strcmp(method,'sum')
    prop_ind=cell2mat(prop_ind.');
    if numel(P_n)==size(prop_ind,1)-numel(unit)
        P_out=sparse(prop_ind(:,1),prop_ind(:,2),[unit; P_n(:)]);
    elseif numel(P_n)*4^(nspins-1)==size(prop_ind,1)
        P_1=repmat(P_in(:,n), 4^(nspins-1),1);
        P_out=sparse(prop_ind(:,1),prop_ind(:,2),P_1(:));
    elseif size(P_n,1)==size(prop_ind,1)
        P_out=spalloc(dim,dim,size(P_in,1)*(dim^(nspins-1)));
        for k=1:nspins
            P_2=1;
            for j=1:nspins
                if j==k
                    P_1=P_in(:,n(k));
                    P_1=sparse(prop_ind(:,1),prop_ind(:,2),[unit; P_1(:)]);
                else
                    P_1=speye(dim^(1/nspins));
                end
                P_2=kron(P_2,P_1);
            end
            P_out=P_out+P_2;
        end
    else
        error('incorrect number of propagator indices')
    end
end

% decided if full or sparse
if nnz(P_out)>dim*dim*dense_matrix, P_out=full(P_out); end

end


function [fdty,fwd_traj,fwd_split,P_n]=qoala_fidelity(param,drift_sys,ctrls,amps,wf,init,targ)

% pull the interaction
interactions=drift_sys.inter_prop;

% pull the state space
state_space=param.space;    basis_set=param.basis;

% dense matrix bound and Hamiltonian dimension
dense_mat=param.sparsity;   dim=param.dim;

% Find out the number of time steps
[nsteps,~]=size(wf);

% scale the waveform with the power levels
wf=amps.*wf;

% number of spins
nspins=param.nspins;

% pull which spin to control
spin_control=param.spin_control;
axes_control=param.ctrl_axes;

% pull the splitting parameters
splitspins_consts=drift_sys.spn_mults;
nsplits=length(splitspins_consts);
q_trott=drift_sys.trotter_number;
ind_spn=drift_sys.ind_spn;  ind_int=drift_sys.ind_int;
num_a=length(ind_spn);      num_b=length(ind_int);

% pull the time increments and Trotterise
dt=repmat(param.pulse_dt/q_trott,[nspins 1]);

% initialise the rotation axes
x=zeros(nsteps*nspins,1); y=zeros(nsteps*nspins,1); z=zeros(nsteps*nspins,1);

% determine the rotation axes for each spin
for k=1:nspins
    k_wf=wf(:,spin_control(k,:));
    if isempty(k_wf), k_wf=zeros(size(k_wf,1),2); end
    x(nsteps*(k-1)+1:nsteps*k)=k_wf(:,1);
    y(nsteps*(k-1)+1:nsteps*k)=k_wf(:,2);
    z(nsteps*(k-1)+1:nsteps*k)=drift_sys.offset(k)*ones(nsteps,1);
end

% expand dt and rotation axis to splittings
dt=reshape((dt*splitspins_consts).',[],1);
x=reshape(repmat(x,[1,nsplits]).',[],1);
y=reshape(repmat(y,[1,nsplits]).',[],1);
z=reshape(repmat(z,[1,nsplits]).',[],1);

% calculate propagator elements
P=rodrigues(state_space,basis_set,ctrls,axes_control,[x, y, z],dt,param.gradops,param.spin_control);

% indices of sparse propagator elements
prop_ind=param.prop_ind;

% define a unit state element
if ismember(basis_set,{'sphten'}), unit=1; else, unit=[]; end
P=P.'; P=[unit*ones(1,size(P,2)); P];

% initialise the effectve propagator
indices=(0:nspins-1)*(nsteps*nsplits);

% initialise propagator cell
P_n=cell(nsplits,nsteps);

% preallocate space for the split trajectories within the time steps
fwd_split=cell(length(ind_spn),nsteps);

% preallocate space for trajectories after splitting, and initialise
fwd_traj=cell(1,nsteps+1);      fwd_traj{1}=init;

if strcmp(state_space,'liouville')

    % propagator forwards
    for n = 1:nsteps

        % split operator method for interaction Hamiltonian
        if size(interactions,1)==1

            % constant interaction
            P_split=(interactions(1,:));
        else

            % time dependent interaction Hamiltonian
            P_split=(interactions(n,:));
        end

        % form and store the propagator at each stage of splitting
        for i=1:nsplits
            P_n{i,n}=ind2propagatorPROP(prop_ind,P,indices+(n-1)*nsplits+i,...
                nspins,unit,dim,'prod',dense_mat);
        end

        % decide sparse or full
        if (nnz(fwd_traj{n}))^2>round(dense_mat*dim*dim)
            fwd_traj{n}=full(fwd_traj{n});
        end

        % first stage is only the coupling propagator exp(couple_1)
        fwd_split{1,n}=(P_split{ind_int(1)}*fwd_traj{n});

        % central stages of exp(couple_bi+1)*exp(spins_ai)
        for i=2:num_a
            fwd_split{i,n}=P_split{ind_int(i)}*(P_n{ind_spn(i-1),n}*fwd_split{i-1,n});
        end

        % final stage stored in the forward trajectory
        fwd_traj{n+1}=P_split{ind_int(num_b)}*(P_n{ind_spn(num_a),n}*fwd_split{num_a,n});

    end

else, error('this version of escalade is only coded for Liouville space');
end

% calculate the fidelity
fdty = real(trace(targ'*fwd_traj{end}))/dim;

end


function [data,fdty,grad]=qoala_gradient(param,drift_sys,ctrls,amps,wf,init,targ,fwd_traj,fwd_split,P_n)

% pull the interaction
interactions=drift_sys.inter_prop;

% pull the state space and basis
state_space=param.space;    basis_set=param.basis;

% dense matrix bound and size of Hamiltonian
dense_mat=param.sparsity;   dim=param.dim;

% Find out the number of time steps and controls
[nsteps,kctrls]=size(wf);

% scale the waveform with the power levels
wf=amps.*wf;

% number of spins
nspins=param.nspins;

% pull which spin to control
spin_control=param.spin_control;
axes_control=param.ctrl_axes;

% pull the splitting parameters
splitspins_consts=drift_sys.spn_mults;
nsplits=length(splitspins_consts);
p_order=drift_sys.split_order;  q_trott=drift_sys.trotter_number;
ind_spn=drift_sys.ind_spn;      ind_int=drift_sys.ind_int;
num_a=length(ind_spn);          num_b=length(ind_int);

% pull the time increments and Trotterise
dt=repmat(param.pulse_dt/q_trott,[nspins 1]);

% initialise the rotation axes
x=zeros(nsteps*nspins,1);
y=zeros(nsteps*nspins,1);
z=zeros(nsteps*nspins,1);

% determine the rotation axes for each spin
for k=1:nspins
    k_wf=wf(:,spin_control(k,:));
    if isempty(k_wf)
        k_wf=zeros(size(k_wf,1),2);
    end
    x(nsteps*(k-1)+1:nsteps*k)=k_wf(:,1);
    y(nsteps*(k-1)+1:nsteps*k)=k_wf(:,2);
    z(nsteps*(k-1)+1:nsteps*k)=drift_sys.offset(k)*ones(nsteps,1);
end

% expand dt and rotation axis to splittings
dt=reshape((dt*splitspins_consts).',[],1);
x=reshape(repmat(x,[1,nsplits]).',[],1);
y=reshape(repmat(y,[1,nsplits]).',[],1);
z=reshape(repmat(z,[1,nsplits]).',[],1);

% Use Euler-Rodrigues for propagator elements
[P,dP]=rodrigues(state_space,basis_set,ctrls,axes_control,[x, y, z],dt,param.gradops,param.spin_control);

for n=1:numel(dP), dP{n}=dP{n}.'; end

% check if in store
if nargin<10

    % indices of sparse propagator elements
    prop_ind=param.prop_ind;

    % define a unit state element
    if ismember(basis_set,{'sphten'}), unit=1; else, unit=[]; end
    P=P.'; P=[unit*ones(1,size(P,2)); P];

    % initialise the effectve propagator
    indices=(0:nspins-1)*(nsteps*nsplits);

    % initialise propagator cell
    P_n=cell(nsplits,nsteps);

    % preallocate space for the split trajectories within the time steps
    fwd_split=cell(length(ind_spn),nsteps);

    % preallocate space for trajectories after splitting, and initialise
    fwd_traj=cell(1,nsteps+1); fwd_traj{1}=init;

    if strcmp(state_space,'liouville')

        % propagator forwards
        for n = 1:nsteps

            % split operator method for interaction Hamiltonian
            if size(interactions,1)==1

                % constant interaction
                P_split=(interactions(1,:));
            else

                % time dependent interaction Hamiltonian
                P_split=(interactions(n,:));
            end

            % form and store the propagator at each stage of splitting
            for i=1:nsplits
                P_n{i,n}=ind2propagatorPROP(prop_ind,P,indices+(n-1)*nsplits+i,...
                    nspins,unit,dim,'prod',dense_mat);
            end

            % decide sparse or full
            if (nnz(fwd_traj{n}))^2>round(dense_mat*dim*dim)
                fwd_traj{n}=full(fwd_traj{n});
            end

            % first stage is only the coupling propagator exp(couple_1)
            fwd_split{1,n}=(P_split{ind_int(1)}*fwd_traj{n});

            % central stages of exp(couple_bi+1)*exp(spins_ai)
            for i=2:num_a
                fwd_split{i,n}=P_split{ind_int(i)}*(P_n{ind_spn(i-1),n}*fwd_split{i-1,n});
            end

            % final stage stored in the forward trajectory
            fwd_traj{n+1}=P_split{ind_int(num_b)}*(P_n{ind_spn(num_a),n}*fwd_split{num_a,n});

        end

    else, error('this version of escalade is only coded for Liouville space');
    end


end

% calculate the fidelity
fdty = real(trace(targ'*fwd_traj{end}))/dim;

% preallocate space for backward split trajectory within the time step
bwd_split=cell(length(ind_spn),1);

% initialise the backward trajectory
bwd_traj=targ;

% preallocate gradient
grad = zeros(nsteps,kctrls);

% preallocate the propagator/D-matrix cell
dP_n=cell(nsplits,1); dP_n(:)={zeros(dim,dim)};

if strcmp(state_space,'liouville')

    for n = nsteps:-1:1

        % split operator method for interaction Hamiltonian
        if size(interactions,1)==1

            % constant interaction
            P_split=(interactions(1,:));
        else

            % time dependent interaction Hamiltonian
            P_split=(interactions(n,:));
        end

        % apply splitting method and store each stage
        if p_order>1

            % first stage of exp(couple_b1)*exp(spins_a1)
            bwd_split{1}=P_n{ind_spn(1),n}'*(P_split{ind_int(1)}'*bwd_traj);

            % central stages of exp(couple_bi)*exp(spins_ai)
            for i=2:num_a
                bwd_split{i}=P_n{ind_spn(i),n}'*(P_split{ind_int(i)}'*bwd_split{i-1});
            end

            % final stage stored in the backward trajectory
            bwd_traj=(P_split{ind_int(num_b)}'*bwd_split{num_a});
        else

            % splitting order < 2 is not symmetric
            bwd_split{1}=P_split{ind_int(num_b)}'*(P_n{ind_spn(num_a),n}'*bwd_traj);
            bwd_traj=bwd_split{1};
        end

        % Allocate gradient column
        grad_col=zeros(1,kctrls);

        for k=1:kctrls

            for i=1:nsplits

                % insert D-matrix into store
                dP_n{i}=reshape(dP{k}(:,(n-1)*nsplits+i),dim,dim);

                % decide on sparsity
                if nnz(dP_n{i})>dense_mat*dim*dim
                    dP_n{i}=full(dP_n{i});
                end

            end

            % (no trotterisation coded here)
            for i=1:num_a
                rho_n=dP_n{ind_spn(i),1}*fwd_split{i,n};
                grad_col(k)=grad_col(k)+amps(k)*real(trace(bwd_split{num_b-i}'*rho_n))/dim;
            end

        end

        % Add to the gradient array
        grad(n,:)=grad_col;

    end

else, error('this version of escalade is only coded for Liouville space');
end

% Return trajectory data
data.trajectory.final=fwd_traj{end};


end