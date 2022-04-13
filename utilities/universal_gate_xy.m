% Adaptive optimal control for state-to-state transfer across an M-particle
% system of coupled, two-level systems.
%
% david.goodwin@chem.ox.ac.uk

function optimpulse=universal_gate_xy(varargin)

options=struct(varargin{:});
if      ~isfield(options,'omega'),          error('must provide ''omega'' as an input')
elseif  ~isfield(options,'target'),         error('must provide ''target'' as an input')
elseif  ~isfield(options,'amplitudes'),     error('must provide ''amplitudes'' as an input')
elseif  ~isfield(options,'duration'),       error('must provide ''duration'' as an input')
elseif  ~isfield(options,'increments'),     error('must provide ''increments'' as an input')
elseif  ~isfield(options,'interaction'),    error('must provide ''interaction'' as an input')
elseif  ~isfield(options,'cartops'),    error('must provide ''cartops'' as an input')
end

% take out the required name-value pair inputs
omega=options.omega;        options=rmfield(options,'omega');
prop_targ=options.target;   options=rmfield(options,'target');
pwr=options.amplitudes;     options=rmfield(options,'amplitudes');
T=options.duration;         options=rmfield(options,'duration');
N=options.increments;       options=rmfield(options,'increments');
H_int=options.interaction;  options=rmfield(options,'interaction');
H_ss=options.cartops;       options=rmfield(options,'cartops');

% Read Optimalisation Parameters
defaultopt = struct('penalties',{{'SNS'}},'init_pulse',2*rand(2*length(pwr),N).'-1,...
			        'tol_x',1e-12,'tol_g',1e-12,'max_iter',100,'prop_cache','carry','prop_zeroed',1e-12);

if (~exist('options','var')) 
    options=defaultopt;
else
    f = fieldnames(defaultopt);
    for i=1:length(f)
        if (~isfield(options,f{i})||(isempty(options.(f{i})))), options.(f{i})=defaultopt.(f{i}); end
    end
end

init_pulse=options.init_pulse;   options=rmfield(options,'init_pulse');

% sparse target propagator
prop_targ=sparse(prop_targ);

% Define control parameters
options.space='liouville';                % state space
options.basis='sphten';                   % basis set
options.operators=H_ss;
options.optimcon_fun=@grape_ugate_qoala;  % Optimal control function
options.drift_sys.interaction=H_int;      % Drift without offsets
options.drift_sys.offset=2*pi*omega(:);   % Drift without coupling
options.prop_targ={prop_targ};            % Target state
options.pulse_dur=T;                      % Pulse duration (s)
options.pulse_nsteps=N;                      % 
options.auxmat_method='taylor';           % taylor much faster than pade
options.nspins=length(omega);             % number of spins in system

% define optimisation parameters
options.method='lbfgs';             % Optimisation method
options.n_grads=25;                 % lbfgs - gradient store number

% parameters for the qoala method
options.pauli_operators=H_ss;
options.ctrl_axes=cell(1,2*size(H_ss,2)/3);
options.operators=cell(1,2*size(H_ss,2)/3);
options.spin_control=zeros(options.nspins,2*size(H_ss,2)/3);
options.pwr_levels=zeros(1,2*size(H_ss,2)/3);
for k=1:size(H_ss,2)/3
    options.ctrl_axes{2*k-1}='x';   options.ctrl_axes{2*k}='y';
    options.pwr_levels(2*k-1)=2*pi*pwr(k);
    options.pwr_levels(2*k)=2*pi*pwr(k);
    for m=1:options.nspins
        if isempty(options.operators{2*k-1}) && ~isempty(H_ss{m,3*k-2})
        options.operators{2*k-1}=H_ss{m,3*k-2};
        elseif ~isempty(options.operators{2*k-1}) && ~isempty(H_ss{m,3*k-2})
        options.operators{2*k-1}=options.operators{2*k-1}+H_ss{m,3*k-2};
        end
        if isempty(options.operators{2*k}) && ~isempty(H_ss{m,3*k-1})
        options.operators{2*k}=H_ss{m,3*k-1};
        elseif ~isempty(options.operators{2*k}) && ~isempty(H_ss{m,3*k-1})
        options.operators{2*k}=options.operators{2*k}+H_ss{m,3*k-1};
        end
        if ~isempty(H_ss{m,3*k-2}), options.spin_control(m,2*k-1)=1; end
        if ~isempty(H_ss{m,3*k-1}), options.spin_control(m,2*k)=1; end
    end
end

% parameters for the qoala method with fidelity check at each iteration

options.drift_sys.splitset=[2 3 4];

param=optimconset(options); % parse
optimpulse=fmaxnewton(param,@optimfun_grape_xy,init_pulse);

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

        % Calculate the penalty (this should be moved to ensemble)
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