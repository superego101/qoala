% Gradient Ascent Pulse Engineering (GRAPE) objective function and 
% gradient. Propagates the system through a user-supplied shaped pulse
% from a given initial propagator and projects the result onto the given 
% final effective propagator. 
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
% The fidelity is returned, along with its gradient and Hessian
% with respect to amplitudes of all control operators at every time step
% of the shaped pulse. Directional propagator derivatives for the gradient
% are calculated with an Auxiliary Matrix method.
%
% Auxiliary Matrix method:          http://dx.doi.org/10.1063/1.4928978
%
% Author:       David L. Goodwin 
% Dev. period:  01/2020--04/2022
% Updated:      13/04/22
%
% Contact:
%   david.goodwin@chem.ox.ac.uk
%   david.goodwin@partner.kit.edu

function [data,fdty,grad]=grape_ugate_auxmat(param,drift_sys,ctrls,amps,wf,init,targ)

% Count the outputs
n_outputs=nargout;

% start the grape timer
grp_time=tic;

% matrix dimension of hamiltonian
dim=param.dim;

% Find out the number of time steps and controls
[nsteps,kctrls]=size(wf);

% scale the waveform with the power levels
wf=amps.*wf;

% dim shorthands
if strcmp(param.space,'hilbert')
    dim_sz=[dim dim]; fwd_traj=init; bwd_trajn=targ';
elseif strcmp(param.space,'liouville') && strcmp(param.ctrl_type,'UR')
    dim_sz=[dim dim]; fwd_traj=init; bwd_trajn=targ;
elseif strcmp(param.space,'liouville')
    dim_sz=dim; fwd_traj=full(init); bwd_trajn=full(targ);
end

if strcmp(param.ctrl_type,'UR')
    nrm=dim;
else
    nrm=1;
end

% Get the timing grid
dt=param.pulse_dt;

% propagation parameters
step_method=param.step_method;
prop_zeroed=param.prop_zeroed;

% gradient methods
if n_outputs>2

    % preallocate gradient
    grad = zeros(nsteps,kctrls);

    % initialise propagator cell
    P_n=cell(1,nsteps);

    dP_rho=cell(1,kctrls); dP_rho(:)={zeros([dim_sz,nsteps])};
end

if strcmp(param.space,'liouville')

    if n_outputs>2

        % propagate forwards
        for n = 1:nsteps

            % Pull the current drift
            L=drift_sys.drift;

            % Add current controls
            for k=1:kctrls
                L=L+wf(n,k)*ctrls{k};
            end

            % Create auxiliary state vector
            aux_vec=[0*fwd_traj; fwd_traj];

            for k=1:kctrls

                % Propagate the auxiliary vector
                aux_mat=[L, ctrls{k}; 0*L, L];
                aux_vector=propagate(param.space,aux_mat,dt(n),step_method,prop_zeroed,aux_vec);

                % pull the derivative
                dP_rho{k}(:,:,n)=aux_vector(1:(end/2),:);
            end

            % Propagate forwards
            fwd_traj=propagate(param.space,L,+dt(n),step_method,prop_zeroed,fwd_traj);

        end

    else

        % propagator forwards
        for n = 1:nsteps

            % Pull the current drift
            L=drift_sys.drift;

            % Add current controls
            for k=1:kctrls
                L=L+wf(n,k)*ctrls{k};
            end

            % Propagate the state forwards
            fwd_traj=propagate(param.space,L,+dt(n),step_method,prop_zeroed,fwd_traj);

        end

    end

    % calculate the fidelity
    if strcmp(param.fidelity,'real')
        fdty = real(trace(targ'*fwd_traj))/nrm;
    elseif strcmp(param.fidelity,'square')
        fdty_overlap = (trace(targ'*fwd_traj))/nrm;
        fdty = fdty_overlap*conj(fdty_overlap);
    end

elseif strcmp(param.space,'hilbert')

    if n_outputs>2

        % propagate forwards
        for n = 1:nsteps

            % Pull the current drift
            L=drift_sys.drift;

            % Add current controls
            for k=1:kctrls
                L=L+wf(n,k)*ctrls{k};
            end

            for k=1:kctrls

                % Propagate the auxiliary vector
                aux_mat=[L, ctrls{k}; 0*L, L];
                prop_dirdiff=propagate(param.space,aux_mat,dt(n),'pade',min([prop_zeroed 1e-14]));
                if k==1, P_n{n}=prop_dirdiff(1:end/2,1:end/2); end

                % pull the derivative
                dP_rho{k}(:,:,n)=(P_n{n}'*prop_dirdiff(1:end/2,1+end/2:end))*fwd_traj;
            end
            
            % Propagate forwards
            fwd_traj=P_n{n}*fwd_traj*P_n{n}';

        end

    else

        % propagator forwards
        for n = 1:nsteps

            % Pull the current drift
            L=drift_sys.drift;

            % Add current controls
            for k=1:kctrls
                L=L+wf(n,k)*ctrls{k};
            end

            % Propagate the state forwards
            fwd_traj=propagate(param.space,L,+dt(n),'taylor',prop_zeroed,fwd_traj);

        end

    end

    % calculate the fidelity
    if strcmp(param.fidelity,'real')
        fdty = real(trace(targ'*fwd_traj))/nrm;
    elseif strcmp(param.fidelity,'square')
        fdty_overlap = (trace(targ'*fwd_traj))/nrm;
        fdty = fdty_overlap*conj(fdty_overlap);
    end

end

% store the final state/effective propagator
data.trajectory.final=fwd_traj;

if strcmp(param.space,'liouville')

    % calculation of gradient elements
    if n_outputs>2

        % loop backwards over timesteps (for efficiency)
        for n = nsteps:-1:1

            % Pull the current drift
            L=drift_sys.drift;

            % calculate the gradient elements
            for k=1:kctrls
                
                if strcmp(param.fidelity,'real')
                    grad(n,k) = amps(k)*real(trace(bwd_trajn'*dP_rho{k}(:,:,n)))/nrm;
                elseif strcmp(param.fidelity,'square')
                    grad_overlap = (trace(bwd_trajn'*dP_rho{k}(:,:,n)))/nrm;
                    grad(n,k)=amps(k)*real(grad_overlap*conj(fdty_overlap)+fdty_overlap*conj(grad_overlap));
                end

                % add the current controls
                L=L+wf(n,k)*ctrls{k};
            end

            % adjoint state
            bwd_trajn=propagate(param.space,L,-dt(n),step_method,prop_zeroed,bwd_trajn);
        end
    end

elseif strcmp(param.space,'hilbert')

    % calculation of gradient elements
    if n_outputs>2

        % loop backwards over timesteps (for efficiency)
        for n = nsteps:-1:1

            % adjoint state
            bwd_trajn=P_n{n}'*bwd_trajn*P_n{n};

            % calculate the gradient elements
            for k=1:kctrls
                grad(n,k) = amps(k)*real(2*trace(bwd_trajn*sparse(dP_rho{k}(:,:,n))));

            end

        end
    end
end

% stop the timer
data.timer=toc(grp_time);

end
