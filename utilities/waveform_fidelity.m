% Function to calculate the fidelity from a pulse waveform.
%
% Author:       David L. Goodwin 
% Dev. period:  01/2020--04/2022
% Updated:      13/04/22
%
% Contact:
%   david.goodwin@chem.ox.ac.uk
%   david.goodwin@partner.kit.edu

function [data,fidelity]=waveform_fidelity(ctrl_sys,drift_sys,controls,...
                                              amps,waveform,rho_init,rho_targ)

% pull the state space
prop_chop=2*eps;
step_method=ctrl_sys.step_method;

% normalisation if UR
if strcmp(ctrl_sys.ctrl_type,'UR'), nrm=ctrl_sys.dim; else, nrm=1; end

% Find out the number of time steps
nsteps=size(waveform,1);
waveform=amps.*waveform;

% Get the timing grid
pulse_dt=ctrl_sys.pulse_dt;
                                         
% Preallocate trajectories
fwd_traj=zeros(([size(rho_init) nsteps+1]));

% Initialise forward trajectories
fwd_traj(:,:,1)=rho_init;

% Run the forward and backward propagation
for n=1:nsteps
    
    % Decide the current drift
    if ~iscell(drift_sys.drift)
        
        % Time-independent drift
        L=drift_sys.drift;
        
    else
        
        % Time-dependent drift
        L=drift_sys.drift{n};
        
    end
    
    % Add current controls to the current drift
    for k=1:numel(controls)
        L=L+waveform(n,k)*controls{k};
    end
    
    % Propagate forwards
    fwd_traj(:,:,n+1)=propagate(ctrl_sys.space,L,+pulse_dt(n),step_method,prop_chop,fwd_traj(:,:,n));
    
end

% calculate the fidelity
if strcmp(ctrl_sys.fidelity,'real')
    fidelity = real(trace(rho_targ'*fwd_traj(:,:,end)))/nrm;
elseif strcmp(ctrl_sys.fidelity,'square')
    fdty_overlap = (trace(rho_targ'*fwd_traj(:,:,end)))/nrm;
    fidelity = fdty_overlap*conj(fdty_overlap);
end

% Return trajectory data
data.trajectory.final=squeeze(fwd_traj(:,:,end));

end
