% This is the function wrapper that updates the parses ctrl_param structure
% and creates an optimal control structure, containing information and 
% parameters needed to perform. Syntax:
%
%   ctrl_sys=optimconox(ctrl_param)
%
% Author:       David L. Goodwin 
% Dev. period:  01/2020--04/2022
% Updated:      13/04/22
%
% Contact:
%   david.goodwin@chem.ox.ac.uk
%   david.goodwin@partner.kit.edu
%
% NOTE: This function is loosely based on the Spinach function optimcon(),
%       and should be referenced accordingly.

function ctrl_sys=optimconset(ctrl_param)

if ischar(ctrl_param) 
    
    % pull help from the helper
    str=helper(ctrl_param);

    % inform the user
    fprintf(1,'%s\n',str{:});
    return
end

% Consistency check
grumble(ctrl_param)

% decide on output to console or to file
if ~isfield(ctrl_param,'output') || ((~isa(ctrl_param.output,'double'))&&(~isa(ctrl_param.output,'char')))
    output=1;
else
    if ctrl_param.output==1
        output=ctrl_param.output;
    else
        output=fopen(ctrl_param.output,'a');
        fprintf(output,'%s\n',' ');
        fprintf(output,'%s\n',datetime('now'));
        fprintf(output,'%s\n',' ');
    end
    ctrl_param=rmfield(ctrl_param,'output');
end
ctrl_sys.output=output;

% Get a unique job identifier
if isfield(ctrl_param,'job_id')
    ctrl_sys.job_id=ctrl_param.job_id;
    ctrl_param=rmfield(ctrl_param,'job_id');
    newid=0;
else
    engine=java.security.MessageDigest.getInstance('MD5');
    engine.update(getByteStreamFromArray([clock feature('getpid')]));
    hashstr=typecast(engine.digest,'uint8');
    ctrl_sys.job_id=sprintf('%.2x',double(hashstr));
    newid=1;
end

% Decide scratch destination
if isfield(ctrl_param,'scratchdir')
    
    % Scratch to a user-specified directory
    ctrl_sys.scratchdir=ctrl_param.scratchdir;
    
    % Parse out
    ctrl_param=rmfield(ctrl_param,'scratchdir');
    
else
    
    % Scratch to the default directory
    S = dbstack();
    rootdir=fileparts(fileparts(which(S(1).file)));
    ctrl_sys.scratchdir=[rootdir filesep 'scratch'];
    
end

% create the sctract directory if it doesn't exist
if ~isfolder(ctrl_sys.scratchdir), mkdir(ctrl_sys.scratchdir); end

% create the job directory in the scratch directory
if ~isfolder([ctrl_sys.scratchdir filesep ctrl_sys.job_id]), mkdir([ctrl_sys.scratchdir filesep ctrl_sys.job_id]); end

% Show the banner
fprintf(output,'%s\n','===========================================');
fprintf(output,'%s\n','=                                         =');
fprintf(output,'%s\n','=               QOALA v1.0                =');
fprintf(output,'%s\n','=                                         =');
fprintf(output,'%s\n','=   M.Foroozandeh, D.L.Goodwin, P.Singh   =');
fprintf(output,'%s\n','=                                         =');
fprintf(output,'%s\n','===========================================');
fprintf(output,'%s\n',' ');
fprintf(output,'%s\n',['scratch location: ' ctrl_sys.scratchdir]);
if newid
    fprintf(output,'%s\n',['new job identifier: ' ctrl_sys.job_id]);
else
    fprintf(output,'%s\n',['inherited job identifier: ' ctrl_sys.job_id]);
end
fprintf(output,'%s\n',' ');

if isfield(ctrl_param,'basis')
    ctrl_sys.basis=ctrl_param.basis;
    ctrl_param=rmfield(ctrl_param,'basis');
else
    error('must provide ctrl_param.basis')
end
fprintf(output,'%s\n',[pad('System basis set',60) ctrl_sys.basis]);

if isfield(ctrl_param,'space') || ~strcmp(ctrl_sys.basis,'sphten')
    ctrl_sys.space=ctrl_param.space;
    ctrl_param=rmfield(ctrl_param,'space');
elseif strcmp(ctrl_sys.basis,'sphten')
    if isfield(ctrl_param,'space') && ~strcmp(ctrl_param.space,'liouville')
        fprintf(output,'%s\n',[pad('ctrl_param.space should be ''liouville'' for ctrl_param.basis=''sphten''',60) 'WARNING']);
        fprintf(output,'%s\n',[pad('... setting ctrl_param.space to ''liouville''',60) 'RESOLVE']);
        ctrl_param=rmfield(ctrl_param,'space');
    elseif isfield(ctrl_param,'space')
        ctrl_param=rmfield(ctrl_param,'space');
    end
    ctrl_sys.space='liouville';
else
    error('must provide ctrl_param.space')
end
fprintf(output,'%s\n',[pad('System state space formalism',60) ctrl_sys.space]);

% 
if isfield(ctrl_param,'mults') && isfield(ctrl_param,'nspins')
    if length(ctrl_param.mults) ~= ctrl_param.nspins
        fprintf(output,'%s\n',[pad('length of ctrl_param.mults does not equal ctrl_param.nspins',60) 'WARNING']);
        fprintf(output,'%s\n',[pad('... setting ctrl_param.nspins to length of ctrl_param.mults',60) 'RESOLVE']);
        ctrl_sys.nspins=length(ctrl_param.mults);
    end
    ctrl_sys.mults=ctrl_param.mults;
    ctrl_param=rmfield(ctrl_param,'mults');
elseif isfield(ctrl_param,'nspins')
    ctrl_sys.nspins=ctrl_param.nspins;
    ctrl_param=rmfield(ctrl_param,'nspins');
    fprintf(output,'%s\n',[pad('ctrl_param.mults not provided',60) 'WARNING']);
    fprintf(output,'%s\n',[pad('... assuming all 2-level systems in ctrl_param.mults',60) 'RESOLVE']);
    ctrl_sys.mults=2*ones(1,ctrl_sys.nspins);
elseif isfield(ctrl_param,'mults')
    ctrl_sys.nspins=length(ctrl_param.mults);
    ctrl_sys.mults=ctrl_param.mults;
    ctrl_param=rmfield(ctrl_param,'mults');
end

% shorthand for this parser
unique_mults=sort(unique(ctrl_sys.mults));

% inform the user
for n=1:length(unique_mults)
    fprintf(output,'%s\n',[pad(['Number of ' int2str(unique_mults(n)) '-level systems'],60) ...
                            int2str(numel(ctrl_sys.mults(ctrl_sys.mults==unique_mults(n))))]);
end

% propagator elements to keep -- very useful with sparse matrices
if ~isfield(ctrl_param,'prop_zeroed'), ctrl_param.prop_zeroed=1e-12; end
ctrl_sys.prop_zeroed=ctrl_param.prop_zeroed;
ctrl_param=rmfield(ctrl_param,'prop_zeroed');

% propagator elements to keep -- very useful with sparse matrices
if ~isfield(ctrl_param,'sparsity'), ctrl_param.sparsity=0.15; end
ctrl_sys.sparsity=ctrl_param.sparsity;
ctrl_param=rmfield(ctrl_param,'sparsity');

% Inform the user
fprintf(output,'%s\n',[pad('Remove propagator elements below',60) ...
                    pad(num2str(ctrl_sys.prop_zeroed,'%0.8g'),20)]);

% Inform the user
fprintf(output,'%s\n',[pad('Propagator sparsity',60) ...
                    pad(num2str(ctrl_sys.sparsity,'%0.8g'),20)]);

% propagator hashing and caching
if ~isfield(ctrl_param,'prop_cache'), ctrl_param.prop_cache='carry'; end
ctrl_sys.prop_cache=ctrl_param.prop_cache;
ctrl_param=rmfield(ctrl_param,'prop_cache');

% Inform the user
if strcmp(ctrl_sys.prop_cache,'store')
    fprintf(output,'%s\n',[pad('Propagator precalculation access',60) pad('stored',20)]);
elseif strcmp(ctrl_sys.prop_cache,'carry')
    fprintf(output,'%s\n',[pad('Propagator precalculation access',60) pad('carry through',20)]);
elseif strcmp(ctrl_sys.prop_cache,'calc')
    fprintf(output,'%s\n',[pad('Propagator precalculation access',60) pad('calculate',20)]);
else
    error(['unknown propagator precaculation method: ' ctrl_sys.prop_cache])
end

% find the optimal control function
ctrl_sys.optimcon_fun=ctrl_param.optimcon_fun;
ctrl_param=rmfield(ctrl_param,'optimcon_fun');

% Inform the user
fprintf(output,'%s\n',[pad('Optimal control function',60) '@' func2str(ctrl_sys.optimcon_fun)]);

% grape algorithm type
if isfield(ctrl_param,'prop_targ')
    
    % universal rotation
    ctrl_sys.ctrl_type='UR';
    
    % check for initial propagator
    if isfield(ctrl_param,'prop_init')

        % grab and parse
        ctrl_sys.initials=ctrl_param.prop_init;   ctrl_param=rmfield(ctrl_param,'prop_init');
    else
        % assume identity as initial propagator
        ctrl_sys.initials=ctrl_param.prop_targ;
        for n=1:numel(ctrl_sys.initials)
            ctrl_sys.initials{n}=speye(size(ctrl_param.prop_targ{n}));
        end
    end
    ctrl_sys.targets=ctrl_param.prop_targ;    ctrl_param=rmfield(ctrl_param,'prop_targ');
    
    
elseif isfield(ctrl_param,'rho_init') && isfield(ctrl_param,'rho_targ')
    
    % assume point2point
    ctrl_sys.ctrl_type='PP';

    % grab the initial and target state pair and parse
    ctrl_sys.initials=ctrl_param.rho_init;  ctrl_param=rmfield(ctrl_param,'rho_init');
    ctrl_sys.targets=ctrl_param.rho_targ;   ctrl_param=rmfield(ctrl_param,'rho_targ');

end

% inform user
if ismember(ctrl_sys.ctrl_type,'PP')
    fprintf(output,'%s\n',[pad('Control algorithm type',60),'point-to-point']);
elseif ismember(ctrl_sys.ctrl_type,'UR')
    fprintf(output,'%s\n',[pad('Control algorithm type',60),'universal rotation']);
end

% Absorb fidelity type
if isfield(ctrl_param,'fidelity')
    ctrl_sys.fidelity=ctrl_param.fidelity;
    ctrl_param=rmfield(ctrl_param,'fidelity');
else
    ctrl_sys.fidelity='real';
end

% Inform the user
switch ctrl_sys.fidelity
    
    case 'real'
        
        % Real part of the overlap
        fprintf(output,'%s\n',[pad('Fidelity measure, range [-1,+1]',60)...
                            pad('Re(<target|rho(T)>)',20)]);
                          
    case 'square'
        
        % Absolute square of the overlap
        fprintf(output,'%s\n',[pad('Fidelity measure, range [0,+1]',60)...
                            pad('|<target|rho(T)>|^2',20)]);
                        
    otherwise
        
        % Complain and bomb out
        error('unknown fidelity functional type.');
                          
end

if isfield(ctrl_param,'fidelity_chk')
    
    % Absorb the specification
    ctrl_sys.fidelity_chk=ctrl_param.fidelity_chk;
    ctrl_param=rmfield(ctrl_param,'fidelity_chk');
    
    % inform the user
    fprintf(output,'%s\n',[pad('fidelity check function',60),'@' func2str(ctrl_sys.fidelity_chk)]);
    
end

% Inform the user
switch ctrl_sys.ctrl_type
    case 'PP'
        fprintf(output,'%s\n',[pad('Initial states per ensemble member',60) ...
            int2str(numel(ctrl_sys.initials))]);
        fprintf(output,'%s\n',[pad('Target states per ensemble member',60) ...
            int2str(numel(ctrl_sys.targets))]);
    case 'UR'
        fprintf(output,'%s\n',[pad('Initial propagators per ensemble member',60) ...
            int2str(numel(ctrl_sys.initials))]);
        fprintf(output,'%s\n',[pad('Target propagators per ensemble member',60) ...
            int2str(numel(ctrl_sys.targets))]);

end

% Absorb control operators
ctrl_sys.operators=ctrl_param.operators;
ctrl_param=rmfield(ctrl_param,'operators');

% pull the dimension of the system
ctrl_sys.dim=length(ctrl_sys.operators{1});
fprintf(output,'%s\n',[pad('Dimension of systems',60) int2str(ctrl_sys.dim)]);

% Inform the user
fprintf(output,'%s\n',[pad('Number of control operators',60) ...
                    int2str(numel(ctrl_sys.operators))]);

% Build the commutation matrix
ctrl_sys.commute=false(numel(ctrl_sys.operators),...
                                  numel(ctrl_sys.operators));
for n=1:numel(ctrl_sys.operators)
    for m=1:numel(ctrl_sys.operators)
        comm=ctrl_sys.operators{n}*ctrl_sys.operators{m}-...
             ctrl_sys.operators{m}*ctrl_sys.operators{n};
        ctrl_sys.commute(n,m)=(norm(comm,1)<1e-6);
    end
end

% Inform the user
fprintf(output,'%s\n',[pad('Number of commuting control pairs',60) int2str((nnz(ctrl_sys.commute)-numel(ctrl_sys.operators))/2)]);
  
if isfield(ctrl_param,'ctrl_axes')

    % Absorb the control axes
    ctrl_sys.ctrl_axes=ctrl_param.ctrl_axes;
    ctrl_param=rmfield(ctrl_param,'ctrl_axes');

    for k=1:numel(ctrl_sys.operators)
        fprintf(output,'%s\n',[pad(['Axis of rotation for control channel ' int2str(k)],60) ...
            [ctrl_sys.ctrl_axes{k} '-axis']]);
    end

else
    ctrl_sys.ctrl_axes={};
end

% Absorb power levels
if iscell(ctrl_param.pwr_levels)
    N=prod(cellfun(@length,ctrl_param.pwr_levels(:)));

    ctrl_sys.pwr_levels=zeros(N,numel(ctrl_sys.operators));
    for n=1:numel(ctrl_sys.operators)
        tmp=1;
        for m=1:numel(ctrl_sys.operators)
            if n==m
                tmp=kron(tmp,(ctrl_param.pwr_levels{m}(:)).');
            else
                tmp=kron(tmp,ones(1,length(ctrl_param.pwr_levels{m})));
            end
        end
        ctrl_sys.pwr_levels(:,n)=tmp;
    end
elseif size(ctrl_param.pwr_levels,2)==1
    ctrl_sys.pwr_levels=ones(1,numel(ctrl_sys.operators)).*ctrl_param.pwr_levels;
elseif size(ctrl_param.pwr_levels,2)==numel(ctrl_sys.operators)
    ctrl_sys.pwr_levels=ctrl_param.pwr_levels;
end

% inform the user
if size(ctrl_sys.pwr_levels,2)>1

    % total simulations
    fprintf(output,'%s\n',[pad('Total number of power level permutations',60) int2str(size(ctrl_sys.pwr_levels,2))]);

    % Inform the user
    for m=1:size(ctrl_sys.pwr_levels,1)
        fprintf(output,'%s\n',[pad(['Unique power levels on control channel ' int2str(m)],60) int2str(length(unique(ctrl_sys.pwr_levels(m,:))))]);
        if length(unique(ctrl_sys.pwr_levels(m,:)))>1
            [val,p]=findprefix(min(ctrl_sys.pwr_levels(m,:)/(2*pi))); fprintf(output,'%s\n',[pad(['Minimum power multiplier for control channel ' int2str(m)],60) [pad(num2str(val,'%7.3f'),7,'right') ' ' p 'Hz']]);
            [val,p]=findprefix(mean(ctrl_sys.pwr_levels(m,:)/(2*pi))); fprintf(output,'%s\n',[pad(['Average power multiplier for control channel ' int2str(m)],60) [pad(num2str(val,'%7.3f'),7,'right') ' ' p 'Hz']]);
            [val,p]=findprefix(max(ctrl_sys.pwr_levels(m,:)/(2*pi))); fprintf(output,'%s\n',[pad(['Maximum power multiplier for control channel ' int2str(m)],60) [pad(num2str(val,'%7.3f'),7,'right') ' ' p 'Hz']]);
        else        
            [val,p]=findprefix(unique(ctrl_sys.pwr_levels(m,:))/(2*pi)); fprintf(output,'%s\n',[pad(['Power multiplier for control channel ' int2str(m)],60) [pad(num2str(val,'%7.3f'),7,'right') ' ' p 'Hz']]);
        end
    end
else
    for m=1:size(ctrl_sys.pwr_levels,1)
        [val,p]=findprefix(ctrl_sys.pwr_levels(m,:)/(2*pi)); fprintf(output,'%s\n',[pad(['Power multiplier for control channel ' int2str(m)],60) [pad(num2str(val,'%7.3f'),7,'right') ' ' p 'Hz']]);
    end
end
ctrl_param=rmfield(ctrl_param,'pwr_levels');


% Process penalties
if isfield(ctrl_param,'penalties')
    
    % Absorb penalties and weights
    ctrl_sys.penalties=ctrl_param.penalties; ctrl_param=rmfield(ctrl_param,'penalties');
    if ~isfield(ctrl_param,'p_weights') && isfield(ctrl_param,'pulse_nsteps')
        ctrl_sys.p_weights=ctrl_param.pulse_nsteps;
    elseif ~isfield(ctrl_param,'p_weights') && isfield(ctrl_param,'pulse_dt')
        ctrl_sys.p_weights=length(ctrl_param.pulse_dt);
    else
        ctrl_sys.p_weights=ctrl_param.p_weights; ctrl_param=rmfield(ctrl_param,'p_weights');
    end
    
else
    
    % Default is no penalty
    ctrl_sys.penalties={'none'};
    ctrl_sys.p_weights=0;
    
end

% Inform the user
for n=1:numel(ctrl_sys.penalties)
    fprintf(output,'%s\n',[pad(['Penalty function ' num2str(n) ', weight ' ...
        num2str(ctrl_sys.p_weights(n),'%.9g')],60)...
        ctrl_sys.penalties{n}]);
end

% Process floor and ceiling bounds
if isfield(ctrl_param,'u_bound')
    
    % Absorb bounds
    ctrl_sys.u_bound=ctrl_param.u_bound; ctrl_param=rmfield(ctrl_param,'u_bound');
    ctrl_sys.l_bound=ctrl_param.l_bound; ctrl_param=rmfield(ctrl_param,'l_bound');
    
else
    
    % Default is nominal power
    ctrl_sys.u_bound=+1;
    ctrl_sys.l_bound=-1;
    
end

% Inform the user
if any(ismember(ctrl_sys.penalties,{'SNS'}))
    fprintf(output,'%s\n',[pad('Maximum SNS penalty ceiling, fraction of power level',60) ...
        num2str(max(ctrl_sys.u_bound(:)),'%+.9g')]);
    fprintf(output,'%s\n',[pad('Minimum SNS penalty floor, fraction of power level',60) ...
        num2str(min(ctrl_sys.l_bound(:)),'%+.9g')]);
end

% Inform the user
if any(ismember(ctrl_sys.penalties,{'SNSA'}))
    fprintf(output,'%s\n',[pad('Maximum SNSA penalty ceiling, fraction of power level',60) ...
        num2str(max(ctrl_sys.u_bound(:)),'%+.9g')]);
end

% zeroed ensemble
if any(ismember(ctrl_sys.penalties,{'ENSNS'}))
    fprintf(output,'%s\n',pad(['Warning: ensemble penalty function ' ctrl_sys.penalties{n} ' not functional'],60));
end

% (non)adiabatic penalty
if any(ismember(ctrl_sys.penalties,{'ADIAB'})) && isfield(ctrl_param,'bandwidth')
    
    ctrl_sys.bandwidth=ctrl_param.bandwidth;
    ctrl_param=rmfield(ctrl_param,'bandwidth');
    
    % save the specification
    fprintf(output,'%s\n',[pad('(Non-)Adiabatic penalty function bandwidth, Hz',60) ...
        num2str(ctrl_sys.bandwidth,'%+.9g')]);
    
elseif any(ismember(ctrl_sys.penalties,{'ADIAB'})) && isfield(ctrl_param,'offsets')
    
    ctrl_sys.bandwidth=range(ctrl_param.offsets);
    
    % save the specification
    fprintf(output,'%s\n',[pad('(Non-)Adiabatic penalty function bandwidth, Hz',60) ...
        num2str(ctrl_sys.bandwidth,'%+.9g')]);
    
elseif any(ismember(ctrl_sys.penalties,{'ADIAB'}))
    
    error('A pulse bandwidth must be provided with the penalty ADIAB')
    
end

if any(ismember(ctrl_sys.penalties,{'ADIAB'})) && isfield(ctrl_sys,'bandwidth')
    
    % create two states for projections
    if isfield(ctrl_param,'adiab_proj') && iscell(ctrl_param.adiab_proj) && numel(ctrl_param.adiab_proj)==3
        ctrl_sys.adiab_proj=ctrl_param.adiab_proj;
        ctrl_param=rmfield(ctrl_param,'adiab_proj');
    else
        error('need to provide x, y, and z states for projection to use ctrl_param.adiab_proj')
    end
    
end

% Absorb the amplitude profile
if isfield(ctrl_param,'amplitudes')
    
    % Take the user-specified values
    ctrl_sys.amplitudes=ctrl_param.amplitudes;
    ctrl_param=rmfield(ctrl_param,'amplitudes');
    fprintf(output,'%s\n',[pad('Optimise only phase of waveform',60) 'yes']);
    
end

% Absorb the amplitude profile
if isfield(ctrl_param,'phases')
    
    % Take the user-specified values
    ctrl_sys.phases=ctrl_param.phases;
    ctrl_param=rmfield(ctrl_param,'phases');
    fprintf(output,'%s\n',[pad('Optimise only amplitude of waveform',60) 'yes']);
    
end

% Absorb pulse timing parameters
if isfield(ctrl_param,'pulse_dt')
    
    % Get the timing grid from the user
    ctrl_sys.pulse_dt=ctrl_param.pulse_dt;
    ctrl_sys.pulse_dur=sum(ctrl_param.pulse_dt);
    ctrl_sys.pulse_nsteps=length(ctrl_param.pulse_dt);
    ctrl_param=rmfield(ctrl_param,'pulse_dt');
    if isfield(ctrl_param,'pulse_dur') && isfield(ctrl_param,'pulse_nsteps')
        fprintf(output,'%s\n',[pad('pulse_dur, pulse_nsteps, and pulse_dt provided',60) 'WARNING']);
        fprintf(output,'%s\n',[pad('... pulse_dt takes precedence',60) 'RESOLVE']);
        ctrl_param=rmfield(ctrl_param,'pulse_dur');
        ctrl_param=rmfield(ctrl_param,'pulse_nsteps');
    elseif isfield(ctrl_param,'pulse_dur')
        fprintf(output,'%s\n',[pad('pulse_dur and pulse_dt provided',60) 'WARNING']);
        fprintf(output,'%s\n',[pad('... pulse_dt takes precedence',60) 'RESOLVE']);
        ctrl_param=rmfield(ctrl_param,'pulse_dur');
    elseif isfield(ctrl_param,'pulse_nsteps')
        fprintf(output,'%s\n',[pad('pulse_nsteps and pulse_dt provided',60) 'WARNING']);
        fprintf(output,'%s\n',[pad('... pulse_dt takes precedence',60) 'RESOLVE']);
        ctrl_param=rmfield(ctrl_param,'pulse_nsteps');
    end
    
else
    
    % Compute the timing grid
    ctrl_sys.pulse_dur=ctrl_param.pulse_dur; 
    ctrl_sys.pulse_nsteps=ctrl_param.pulse_nsteps; 
    ctrl_sys.pulse_dt=ones(ctrl_param.pulse_nsteps,1)*(ctrl_sys.pulse_dur/...
                                                               ctrl_sys.pulse_nsteps);
    ctrl_param=rmfield(ctrl_param,'pulse_dur'); ctrl_param=rmfield(ctrl_param,'pulse_nsteps');
    
end

if length(unique(ctrl_sys.pulse_dt))==1
    fprintf(output,'%s\n',[pad('Waveform slice duration grid',60) 'uniform']);
    [val,p]=findprefix(ctrl_sys.pulse_dt(1));
    fprintf(output,'%s\n',[pad('Control sequence slice duration',60) [pad(num2str(val,'%7.3f'),7,'right') ' ' p 's']]);
else
    fprintf(output,'%s\n',[pad('Waveform slice duration grid',60) 'non-uniform']);
    [val,p]=findprefix(mean(ctrl_sys.pulse_dt));
    fprintf(output,'%s\n',[pad('Average control sequence slice duration',60) [pad(num2str(val,'%7.3f'),7,'right') ' ' p 's']]);

    fprintf(output,'%s\n',[pad('non-uniform slice duration grid not fully coded for operator splitting',60) 'WARNING']);
end

% Inform the user
fprintf(output,'%s\n',[pad('Number of slices in the control sequence',60) ...
                    int2str(ctrl_sys.pulse_nsteps)]);
[val,p]=findprefix(ctrl_sys.pulse_dur);
fprintf(output,'%s\n',[pad('Total duration of the control sequence',60) [pad(num2str(val,'%7.3f'),7,'right') ' ' p 's']]);

% drift systems
if isfield(ctrl_param,'drift_sys') && iscell(ctrl_param.drift_sys)

    % Inform the user
    fprintf(output,'%s\n',[pad('Ensemble size',60) ...
                    int2str(numel(ctrl_sys.drift_sys))]);
elseif isfield(ctrl_param,'drift_sys')
    % ensure conistencey, and covert to a cell
    drift_sys=ctrl_param.drift_sys;
    ctrl_param=rmfield(ctrl_param,'drift_sys');
    ctrl_param.drift_sys{1}=drift_sys;
else
    ctrl_param.drift_sys=cell(1,1);
end

% ensure cell arrays for all drift system parameters
for n=1:numel(ctrl_param.drift_sys)
    if numel(ctrl_param.drift_sys)>1
        ens_ind=['{' int2str(n) '}'];
        sys_str=['Ensemble member ' ens_ind ': '];
    else
        ens_ind='';
        sys_str='';
    end
    drift_sys_in=ctrl_param.drift_sys{n};

    % inform the user, and set correct parameters for escalade or auxmat
    if contains(func2str(ctrl_sys.optimcon_fun),{'escalade','qoala'}) &&...
        (isfield(ctrl_sys,'fidelity_chk') && contains(func2str(ctrl_sys.fidelity_chk),'waveform_fidelity'))
        if ~isfield(drift_sys_in,'offset')
            error([ sys_str ' offset field required for @' func2str(ctrl_sys.optimcon_fun)]);
        end
        drift_sys_out.offset=drift_sys_in.offset;
        if isfield(drift_sys_in,'interaction')
            drift_sys_out.interaction=drift_sys_in.interaction;
            fprintf(output,'%s\n',[pad([sys_str 'offset and interaction will be used for @' func2str(ctrl_sys.optimcon_fun)],60)]);
        else
            fprintf(output,'%s\n',[pad([sys_str 'offset will be used for @' func2str(ctrl_sys.optimcon_fun)],60)]);
        end
        if isfield(drift_sys_in,'interaction') && isfield(drift_sys_in,'singlespin')
            drift_sys_out.drift=drift_sys_in.interaction+drift_sys_in.singlespin;
            fprintf(output,'%s\n',[pad([sys_str 'singlespin and interaction will be used for @' func2str(ctrl_sys.fidelity_chk)],60)]);
        elseif isfield(drift_sys_in,'interaction')
            drift_sys_out.drift=drift_sys_in.interaction;
            fprintf(output,'%s\n',[pad([sys_str 'interaction will be used for @' func2str(ctrl_sys.fidelity_chk)],60)]);
        elseif isfield(drift_sys_in,'singlespin')
            drift_sys_out.drift=drift_sys_in.singlespin;
            fprintf(output,'%s\n',[pad([sys_str 'singlespin will be used for @' func2str(ctrl_sys.fidelity_chk)],60)]);
        elseif nnz(drift_sys_out.offset)==0
            drift_sys_out.drift=sparse(zeros(ctrl_sys.dim,ctrl_sys.dim));
            fprintf(output,'%s\n',[pad([sys_str 'interaction and/or singlespin required for @' func2str(ctrl_sys.fidelity_chk)],60) 'WARNING']);
            fprintf(output,'%s\n',[pad([sys_str '... assuming interaction + singlespin = 0'],60) 'RESOLVE']);
        else
            error(['non-zeros offsets detected in drift_sys.offsets, but not supplied in matrix form for @' func2str(ctrl_sys.fidelity_chk)])
        end

    elseif contains(func2str(ctrl_sys.optimcon_fun),{'escalade','qoala'})
        if ~isfield(drift_sys_in,'offset')
            error([ sys_str ' offset field required for @' func2str(ctrl_sys.optimcon_fun)]);
        end
        drift_sys_out.offset=drift_sys_in.offset;
        if isfield(drift_sys_in,'interaction')
            drift_sys_out.interaction=drift_sys_in.interaction;
            fprintf(output,'%s\n',[pad([sys_str 'offset and interaction will be used for @' func2str(ctrl_sys.optimcon_fun)],60)]);
        else
            fprintf(output,'%s\n',[pad([sys_str 'offset will be used for @' func2str(ctrl_sys.optimcon_fun)],60)]);
        end

    else
        if isfield(drift_sys_in,'interaction') && isfield(drift_sys_in,'singlespin')
            drift_sys_out.drift=drift_sys_in.interaction+drift_sys_in.singlespin;
            fprintf(output,'%s\n',[pad([sys_str 'singlespin and interaction will be used for @' func2str(ctrl_sys.optimcon_fun)],60)]);
        elseif isfield(drift_sys_in,'interaction')
            drift_sys_out.drift=drift_sys_in.interaction;
            fprintf(output,'%s\n',[pad([sys_str 'interaction will be used for @' func2str(ctrl_sys.optimcon_fun)],60)]);
        elseif isfield(drift_sys_in,'singlespin')
            drift_sys_out.drift=drift_sys_in.singlespin;
            fprintf(output,'%s\n',[pad([sys_str 'singlespin will be used for @' func2str(ctrl_sys.optimcon_fun)],60)]);
        else
            drift_sys_out.drift=sparse(zeros(ctrl_sys.dim,ctrl_sys.dim));
            fprintf(output,'%s\n',[pad([sys_str 'interaction and/or singlespin required for @' func2str(ctrl_sys.optimcon_fun)],60) 'WARNING']);
            fprintf(output,'%s\n',[pad([sys_str '... assuming interaction + singlespin = 0'],60) 'RESOLVE']);
        end

    end

    % check for time dependence and for zeros in interaction
    if isfield(drift_sys_out,'drift')
        if iscell(drift_sys_out.drift)
            if (all(cellfun(@nnz,drift_sys_out.drift(:))==0))
                fprintf(output,'%s\n',[pad([sys_str 'all drift(t) = 0 detected'],60) 'WARNING']);
                drift_sys_out.drift=0;
                fprintf(output,'%s\n',[pad([sys_str '... drift time-dependence removed'],60) 'RESOLVE']);
                fprintf(output,'%s\n',[pad([sys_str '2-norm of drift'],60) [pad(num2str(0,'%7.3f'),7,'right') ' Hz']]);
            else
                fprintf(output,'%s\n',[pad([sys_str 'Time-dependent drift'],60) 'yes']);
                [val,p]=findprefix((mean(cellfun(@norm,cellfun(@full,(drift_sys_out.drift(:)),'UniformOutput',0))))/(2*pi));
                fprintf(output,'%s\n',[pad([sys_str 'average 2-norm of drifts'],60) [pad(num2str(val,'%7.3f'),7,'right') ' ' p 'Hz']]);
            end
        else
            if nnz(drift_sys_out.drift)==0
                drift_sys_out.drift=0;
            end
            fprintf(output,'%s\n',[pad([sys_str 'Time-dependent drift'],60) 'no']);
            [val,p]=findprefix(norm(full(drift_sys_out.drift))/(2*pi));
            fprintf(output,'%s\n',[pad([sys_str '2-norm of drift'],60) [pad(num2str(val,'%7.3f'),7,'right') ' ' p 'Hz']]);
        end
    end
    if isfield(drift_sys_out,'offset')
        if iscell(drift_sys_out.offset)
            if (all(cellfun(@nnz,drift_sys_out.offset(:))==0))
                fprintf(output,'%s\n',[pad([sys_str 'all offset(t) = 0 detected'],60) 'WARNING']);
                ctrl_sys_chk.offset=ctrl_sys_chk.offset{1};
                fprintf(output,'%s\n',[pad([sys_str '... offset time-dependence removed'],60) 'RESOLVE']);
                fprintf(output,'%s\n',[pad([sys_str '2-norm of offset'],60) [pad(num2str(0,'%7.3f'),7,'right') ' Hz']]);
            else
                fprintf(output,'%s\n',[pad([sys_str 'Time-dependent offset'],60) 'yes']);
                [val,p]=findprefix((mean(cellfun(@norm,cellfun(@full,(drift_sys_out.offset(:)),'UniformOutput',0))))/(2*pi));
                fprintf(output,'%s\n',[pad([sys_str 'average 2-norm of offsets'],60) [pad(num2str(val,'%7.3f'),7,'right') ' ' p 'Hz']]);
            end
        else
            fprintf(output,'%s\n',[pad([sys_str 'Time-dependent offset'],60) 'no']);            
            [val,p]=findprefix(norm(full(drift_sys_out.offset))/(2*pi));
            fprintf(output,'%s\n',[pad([sys_str '2-norm of offset'],60) [pad(num2str(val,'%7.3f'),7,'right') ' ' p 'Hz']]);
        end
    end
    if isfield(drift_sys_out,'singlespin')
        if iscell(drift_sys_out.singlespin)
            if (all(cellfun(@nnz,drift_sys_out.singlespin(:))==0))
                fprintf(output,'%s\n',[pad([sys_str 'all singlespin(t) = 0 detected'],60) 'WARNING']);
                drift_sys_out.singlespin=0;
                fprintf(output,'%s\n',[pad([sys_str '... singlespin time-dependence removed'],60) 'RESOLVE']);
                fprintf(output,'%s\n',[pad([sys_str '2-norm of singlespin'],60) [pad(num2str(0,'%7.3f'),7,'right') ' Hz']]);
            else
                fprintf(output,'%s\n',[pad([sys_str 'Time-dependent singlespin'],60) 'yes']);
                [val,p]=findprefix((mean(cellfun(@norm,cellfun(@full,(drift_sys_out.spinglespin(:)),'UniformOutput',0))))/(2*pi));
                fprintf(output,'%s\n',[pad([sys_str 'average 2-norm of singlespins'],60) [pad(num2str(val,'%7.3f'),7,'right') ' ' p 'Hz']]);
            end
        else
            if nnz(drift_sys_out.singlespin)==0
                drift_sys_out.singlespin=0;
            end
            fprintf(output,'%s\n',[pad([sys_str 'Time-dependent singlespin'],60) 'no']);            
            [val,p]=findprefix(norm(full(drift_sys_out.singlespin))/(2*pi));
            fprintf(output,'%s\n',[pad([sys_str '2-norm of singlespin'],60) [pad(num2str(val,'%7.3f'),7,'right') ' ' p 'Hz']]);
        end
    end
    if isfield(drift_sys_out,'interaction')
        if iscell(drift_sys_out.interaction)
            if (all(cellfun(@nnz,drift_sys_out.interaction(:))==0))
                fprintf(output,'%s\n',[pad([sys_str 'all interaction(t) = 0 detected'],60) 'WARNING']);
                drift_sys_out=rmfield(drift_sys_out,'interaction');
                fprintf(output,'%s\n',[pad([sys_str '... interaction removed from ctrl_sys.drift_sys' ens_ind],60) 'RESOLVE']);
            else
                fprintf(output,'%s\n',[pad([sys_str 'Time-dependent interaction'],60) 'yes']);
                [val,p]=findprefix((mean(cellfun(@norm,cellfun(@full,(drift_sys_out.interaction(:)),'UniformOutput',0))))/(2*pi));
                fprintf(output,'%s\n',[pad([sys_str 'average 2-norm of interactions'],60) [pad(num2str(val,'%7.3f'),7,'right') ' ' p 'Hz']]);
            end            
        else
            if nnz(drift_sys_out.interaction)==0
                fprintf(output,'%s\n',[pad([sys_str 'interaction = 0 detected'],60) 'WARNING']);
                fprintf(output,'%s\n',[pad([sys_str '... interaction removed from ctrl_sys.drift_sys' ens_ind],60) 'RESOLVE']);
                drift_sys_out=rmfield(drift_sys_out,'interaction');
            else
                fprintf(output,'%s\n',[pad([sys_str 'Time-dependent interaction'],60) 'no']);            
                [val,p]=findprefix(norm(full(drift_sys_out.interaction))/(2*pi));
                fprintf(output,'%s\n',[pad([sys_str '2-norm of interaction'],60) [pad(num2str(val,'%7.3f'),7,'right') ' ' p 'Hz']]);
            end
        end
    end
    if ~isfield(drift_sys_out,'interaction')
        if isfield(drift_sys_in,'split_method') ||...
                isfield(drift_sys_in,'split_order') ||...
                isfield(drift_sys_in,'trotter_number') ||...
                isfield(drift_sys_in,'trotterset') ||...
                isfield(drift_sys_in,'adaptset') ||...
                isfield(drift_sys_in,'splitset') ||...
                isfield(drift_sys_in,'adapt_tols') ||...
                isfield(drift_sys_in,'adapt_method') ||...
                isfield(drift_sys_in,'adapt_scale') ||...
                isfield(drift_sys_in,'adapt_miniter') ||...
                isfield(drift_sys_in,'adapt_interval')
            fprintf(output,'%s\n',[pad([sys_str 'No interaction present'],60) 'WARNING']);
        end
        if isfield(drift_sys_in,'split_method'),   fprintf(output,'%s\n',[pad([sys_str '... split_method removed from ctrl_sys.drift_sys' ens_ind],60) 'RESOLVE']); end
        if isfield(drift_sys_in,'split_order'),    fprintf(output,'%s\n',[pad([sys_str '... split_order removed from ctrl_sys.drift_sys' ens_ind],60) 'RESOLVE']); end
        if isfield(drift_sys_in,'trotter_number'), fprintf(output,'%s\n',[pad([sys_str '... trotter_number removed from ctrl_sys.drift_sys' ens_ind],60) 'RESOLVE']); end
        if isfield(drift_sys_in,'trotterset'), fprintf(output,'%s\n',[pad([sys_str '... trotterset removed from ctrl_sys.drift_sys' ens_ind],60) 'RESOLVE']); end
        if isfield(drift_sys_in,'adaptset'),   fprintf(output,'%s\n',[pad([sys_str '... adaptset removed from ctrl_sys.drift_sys' ens_ind],60) 'RESOLVE']); end
        if isfield(drift_sys_in,'splitset'),   fprintf(output,'%s\n',[pad([sys_str '... splitset removed from ctrl_sys.drift_sys' ens_ind],60) 'RESOLVE']); end
        if isfield(drift_sys_in,'adapt_tols'),     fprintf(output,'%s\n',[pad([sys_str '... adapt_tols removed from ctrl_sys.drift_sys' ens_ind],60) 'RESOLVE']); end
        if isfield(drift_sys_in,'adapt_method'),   fprintf(output,'%s\n',[pad([sys_str '... adapt_method removed from ctrl_sys.drift_sys' ens_ind],60) 'RESOLVE']); end
        if isfield(drift_sys_in,'adapt_scale'),    fprintf(output,'%s\n',[pad([sys_str '... adapt_scale removed from ctrl_sys.drift_sys' ens_ind],60) 'RESOLVE']); end
        if isfield(drift_sys_in,'adapt_miniter'),  fprintf(output,'%s\n',[pad([sys_str '... adapt_miniter removed from ctrl_sys.drift_sys' ens_ind],60) 'RESOLVE']); end
        if isfield(drift_sys_in,'adapt_interval'),  fprintf(output,'%s\n',[pad([sys_str '... adapt_interval removed from ctrl_sys.drift_sys' ens_ind],60) 'RESOLVE']); end
        if ~isfield(ctrl_param,'gradops'), ctrl_param.gradops='sparse'; end

    else
        % detect the splitting method
        if isfield(drift_sys_in,'split_method')
            if      ismember(drift_sys_in.split_method,{'uncoupled','none'}),      split_order=0;
            elseif  ismember(drift_sys_in.split_method,{'Lie','Lie-1st'}),         split_order=1;
            elseif  ismember(drift_sys_in.split_method,{'Strang','Strang-2nd'}),   split_order=2;
            elseif  ismember(drift_sys_in.split_method,{'Blanes-3rd'}),            split_order=3;
            elseif  ismember(drift_sys_in.split_method,{'Blanes','Blanes-4th'}),   split_order=4;
            elseif  ismember(drift_sys_in.split_method,{'Omelyan','Omelyan-6th'}), split_order=6;
            else,   error([sys_str 'Unknown splitting method: ' drift_sys_in.split_method])
            end
            drift_sys_in=rmfield(drift_sys_in,'split_method');
            if isfield(drift_sys_in,'split_order')
                fprintf(output,'%s\n',[pad([sys_str 'split_order and split_method supplied'],60) 'WARNING']);
                fprintf(output,'%s\n',[pad([sys_str '... split_method takes precedence'],60) 'RESOLVE']);
                drift_sys_in=rmfield(drift_sys_in,'split_order');
            end
            if isfield(drift_sys_in,'splitset')
                splitset=drift_sys_in.splitset;
                drift_sys_in=rmfield(drift_sys_in,'splitset');
            end
        elseif isfield(drift_sys_in,'split_order')
            split_order=drift_sys_in.split_order;
            drift_sys_in=rmfield(drift_sys_in,'split_order');
            if isfield(drift_sys_in,'splitset')
                splitset=drift_sys_in.splitset;
                drift_sys_in=rmfield(drift_sys_in,'splitset');
            end
        elseif isfield(drift_sys_in,'splitset')
            splitset=drift_sys_in.splitset;
            drift_sys_in=rmfield(drift_sys_in,'splitset');
            split_order=min(splitset);
        else
            split_order=2;
        end

        % give the splitting mehtod an integer order
        if      split_order==0, split_method='uncoupled';
        elseif  split_order==1, split_method='Lie-1st';
        elseif  split_order==2, split_method='Strang-2nd';
        elseif  split_order==3, split_method='Blanes-3rd';
        elseif  split_order==4, split_method='Blanes-4th';
        elseif  split_order==6, split_method='Omelyan-6th';
        else,   error('unknown split order')
        end


        % decided if Trotterised splittting
        if isfield(drift_sys_in,'trotter_number')
            trotter_number=drift_sys_in.trotter_number;
            drift_sys_in=rmfield(drift_sys_in,'trotter_number');
            if isfield(drift_sys_in,'trotterset')
                trotterset=drift_sys_in.trotterset;
                drift_sys_in=rmfield(drift_sys_in,'trotterset');
            end
        elseif isfield(drift_sys_in,'trotterset')
            trotterset=drift_sys_in.trotterset;
            trotter_number=min(drift_sys_in.trotterset);
        else
            trotter_number=1;
        end

        if exist('splitset','var') && exist('trotterset','var')

            adapt=true;
            if isfield(drift_sys_in,'adapt_tols')
                adapt_tols=drift_sys_in.adapt_tols;
                drift_sys_in=rmfield(drift_sys_in,'adapt_tols');
            else
                adapt_tols=[0.5 1e-8];
            end
            if isfield(drift_sys_in,'adaptset')
                fprintf(output,'%s\n',[pad([sys_str 'adaptset, trotterset and splitset provided'],60) 'WARNING']);
                fprintf(output,'%s\n',[pad([sys_str '... trotterset and splitset take precedence'],60) 'RESOLVE']);
                drift_sys_in=rmfield(drift_sys_in,'adaptset');
            end
            adaptset=[kron(ones(1,length(splitset)),trotterset);...
                kron(splitset,ones(1,length(trotterset)))];

            clear splitset trotterset

        elseif exist('trotterset','var')

            adapt=true;
            if isfield(drift_sys_in,'adapt_tols')
                adapt_tols=drift_sys_in.adapt_tols;
                drift_sys_in=rmfield(drift_sys_in,'adapt_tols');
            else
                adapt_tols=[0.5 1e-8];
            end
            if isfield(drift_sys_in,'adaptset')
                fprintf(output,'%s\n',[pad([sys_str 'adaptset and trotterset provided'],60) 'WARNING']);
                fprintf(output,'%s\n',[pad([sys_str '... trotterset takes precedence'],60) 'RESOLVE']);
                drift_sys_in=rmfield(drift_sys_in,'adaptset');
            end
            adaptset=[trotterset; split_order*ones(1,length(trotterset))];

            clear('trotterset');

        elseif exist('splitset','var')

            adapt=true;
            if isfield(drift_sys_in,'adapt_tols')
                adapt_tols=drift_sys_in.adapt_tols;
                drift_sys_in=rmfield(drift_sys_in,'adapt_tols');
            else
                adapt_tols=[0.5 1e-8];
            end
            if isfield(drift_sys_in,'adaptset')
                fprintf(output,'%s\n',[pad([sys_str 'adaptset and splitset provided'],60) 'WARNING']);
                fprintf(output,'%s\n',[pad([sys_str '... splitset takes precedence'],60) 'RESOLVE']);
                drift_sys_in=rmfield(drift_sys_in,'adaptset');
            end
            adaptset=[trotter_number*ones(1,length(splitset)); splitset];

            clear('splitset');

        elseif isfield(drift_sys_in,'adaptset')

            adapt=true;
            adaptset=drift_sys_in.adaptset;
            drift_sys_in=rmfield(drift_sys_in,'adaptset');
            if isfield(drift_sys_in,'adapt_tols')
                adapt_tols=drift_sys_in.adapt_tols;
                drift_sys_in=rmfield(drift_sys_in,'adapt_tols');
            else
                adapt_tols=[0.5 1e-8];
            end

        else
            adapt=false;
            adaptset=[trotter_number; split_order];
            adapt_tols=[inf -inf];
        end

        if isfield(drift_sys_in,'adapt_method')

            adapt_method=drift_sys_in.adapt_method;
            drift_sys_in=rmfield(drift_sys_in,'adapt_method');
        else
            adapt_method='gain';
        end

        if isfield(drift_sys_in,'adapt_scale')

            adapt_scale=drift_sys_in.adapt_scale;
            drift_sys_in=rmfield(drift_sys_in,'adapt_scale');
        else
            adapt_scale='infidelity';
        end

        if isfield(drift_sys_in,'adapt_minit')

            adapt_minit=drift_sys_in.adapt_minit;
            adapt_counter=inf;
            drift_sys_in=rmfield(drift_sys_in,'adapt_minit');
        else
            adapt_minit=5;
            adapt_counter=inf;
        end

        if isfield(drift_sys_in,'adapt_interval')

            adapt_interval=drift_sys_in.adapt_interval;
            drift_sys_in=rmfield(drift_sys_in,'adapt_interval');
        else
            adapt_interval=1;
        end

        % get the splits
        dt_in=ctrl_sys.pulse_dur/(ctrl_sys.pulse_nsteps); % note that this is for uniform time slicing
        splitting_param.prop_method='taylor';
        splitting_param.space=ctrl_sys.space;
        splitting_param.prop_zeroed=ctrl_sys.prop_zeroed;
        splitting_param.scratchdir=ctrl_sys.scratchdir;
        splitting_param.job_id=ctrl_sys.job_id;
        splitting_param.prop_cache=ctrl_sys.prop_cache;

        switch ctrl_sys.prop_cache

            case {'store'}

                for k=1:size(adaptset,2)

                    if adaptset(2,k)==split_order && adaptset(1,k)==trotter_number

                        [spn_mults,ind_spn,ind_int,P]=operator_splittings(splitting_param,adaptset(2,k),adaptset(1,k),dt_in,drift_sys_out.interaction);
                        drift_sys_out.spn_mults=spn_mults;
                        drift_sys_out.ind_spn=ind_spn;
                        drift_sys_out.ind_int=ind_int;
                        drift_sys_out.inter_prop=P;

                    else, operator_splittings(splitting_param,adaptset(2,k),adaptset(1,k),dt_in,drift_sys_out.interaction);
                    end

                end

            case {'carry'}

                splitting_param.prop_cache='calc';

                for k=1:size(adaptset,2)

                    [spn_mults,ind_spn,ind_int,P]=operator_splittings(splitting_param,adaptset(2,k),adaptset(1,k),dt_in,drift_sys_out.interaction);
                    drift_sys_out.spn_mults{k,1}=spn_mults;
                    drift_sys_out.ind_spn{k,1}=ind_spn;
                    drift_sys_out.ind_int{k,1}=ind_int;
                    drift_sys_out.inter_prop(k,1:numel(P))=P;

                end

                if size(adaptset,2)==1

                    drift_sys_out.spn_mults=drift_sys_out.spn_mults{1};
                    drift_sys_out.ind_spn=drift_sys_out.ind_spn{1};
                    drift_sys_out.ind_int=drift_sys_out.ind_int{1};

                end

            case {'calc'}

                [spn_mults,ind_spn,ind_int,P]=operator_splittings(splitting_param,split_order,trotter_number,dt_in,drift_sys_out.interaction);
                drift_sys_out.spn_mults=spn_mults;
                drift_sys_out.ind_spn=ind_spn;
                drift_sys_out.ind_int=ind_int;
                drift_sys_out.inter_prop=P;
        end

        if size(adaptset,2)>1 && adapt
            if unique(adaptset(2,:))>1
                fprintf(output,'%s\n',[pad([sys_str 'Initial operator splitting method'],60),split_method]);
                fprintf(output,'%s\n',[pad([sys_str 'Initial operator splitting order'],60),int2str(split_order)]);
                fprintf(output,'%s\n',[pad([sys_str 'Minimum operator splitting order'],60),int2str(min(adaptset(2,:)))]);
                fprintf(output,'%s\n',[pad([sys_str 'Maximum operator splitting order'],60),int2str(max(adaptset(2,:)))]);
            else
                fprintf(output,'%s\n',[pad([sys_str 'Operator splitting method'],60),split_method]);
                fprintf(output,'%s\n',[pad([sys_str 'Operator splitting order'],60),int2str(split_order)]);
            end
            if unique(adaptset(1,:))>1
                fprintf(output,'%s\n',[pad([sys_str 'Initial Trotter number'],60),int2str(trotter_number)]);
                fprintf(output,'%s\n',[pad([sys_str 'Minimum Trotter number'],60),int2str(min(adaptset(1,:)))]);
                fprintf(output,'%s\n',[pad([sys_str 'Maximum Trotter number'],60),int2str(max(adaptset(1,:)))]);
            else
                fprintf(output,'%s\n',[pad([sys_str 'Trotter number'],60),int2str(trotter_number)]);
            end
        else
            fprintf(output,'%s\n',[pad([sys_str 'Operator splitting method'],60),split_method]);
            fprintf(output,'%s\n',[pad([sys_str 'Operator splitting order'],60),int2str(split_order)]);
            fprintf(output,'%s\n',[pad([sys_str 'Trotter number'],60),int2str(trotter_number)]);
        end
        
        if size(adaptset,2)>1 && adapt
            fprintf(output,'%s\n',[pad([sys_str 'Minimum fidelity error threshold'],60),num2str(adapt_tols(2),'%+.9g')]);
            fprintf(output,'%s\n',[pad([sys_str 'Maximum fidelity error threshold'],60),num2str(adapt_tols(1),'%+.9g')]);
        elseif size(adaptset,2)==1 && adapt
            fprintf(output,'%s\n',[pad([sys_str 'Fidelity error threshold'],60),num2str(adapt_tols(1),'%+.9g')]);
        end
        if adapt_minit>0 && adapt
            fprintf(output,'%s\n',[pad([sys_str 'Minimum iterations between adaptivity'],60),int2str(adapt_minit)]);
        end
        if adapt_interval>0 && isfinite(adapt_interval) && adapt
            fprintf(output,'%s\n',[pad([sys_str 'Maximum jumps within adaptivity interations'],60),int2str(adapt_interval)]);
        end

        drift_sys_out.trotter_number=trotter_number;
        drift_sys_out.split_order=split_order;
        drift_sys_out.adapt=adapt;
        if adapt
            drift_sys_out.adaptset=adaptset;
            drift_sys_out.adapt_tols=adapt_tols;
            drift_sys_out.adapt_method=adapt_method;
            drift_sys_out.adapt_scale=adapt_scale;
            drift_sys_out.adapt_minit=adapt_minit;
            drift_sys_out.adapt_interval=adapt_interval;
            drift_sys_out.adapt_counter=adapt_counter;
        end
        if split_order~=0 && isfield(ctrl_param,'gradops') && strcmp(ctrl_param.gradops,'sparse')
            fprintf(output,'%s\n',[pad([sys_str 'ctrl_param.gradops=''sparse'' and split_order>0'],60) 'WARNING']);
            fprintf(output,'%s\n',[pad([sys_str '... changing to ctrl_param.gradops=''full'''],60) 'RESOLVE']);
            ctrl_param.gradops='full';
        end

        clear('split_method','split_order','trotter_number',...
              'adapt','adaptset','adapt_tols','adapt_method',...
              'adapt_scale','adapt_minit','adapt_interval','adapt_counter');
    end
    % Catch unparsed fields
    unparsed=fieldnames(drift_sys_in);
    if ~isempty(unparsed)
        for m=1:numel(unparsed)
            fprintf(output,'%s\n',[pad([sys_str 'Unrecognised ctrl_param.drift_sys option - ' unparsed{m}],60) 'WARNING']);
        end
    end
    ctrl_param.drift_sys{n}=drift_sys_out;
    clear('drift_sys_in','drift_sys_out');
end

% absorb the drift ensemble
ctrl_sys.drift_sys=ctrl_param.drift_sys;
ctrl_param=rmfield(ctrl_param,'drift_sys');

% Process cavity effect
if isfield(ctrl_param,'cavity_decay_rate')
    
    % Absorb cavity decay rate
    ctrl_sys.cavity_decay_rate=ctrl_param.cavity_decay_rate;
    ctrl_param=rmfield(ctrl_param,'cavity_decay_rate');
    
    % Inform the user
    fprintf(1,'%s\n',[pad('Cavity decay rate, Hertz',60) ...
        num2str(ctrl_sys.cavity_decay_rate,'%.9g')]);
    
    if ctrl_sys.cavity_decay_rate~=0
        
        % Absorb cavity decay rate
        ctrl_sys.cavity_function=ctrl_param.cavity_function;
        ctrl_param=rmfield(ctrl_param,'cavity_function');
        
        % Inform the user
        fprintf(output,'%s\n',[pad('Pulse cavity effect function',60) ...
            '@' func2str(ctrl_sys.cavity_function)]);
        
        % interpolation points for cavity effect
        ctrl_sys.cavity_n_interp=ctrl_param.cavity_n_interp;
        ctrl_param=rmfield(ctrl_param,'cavity_n_interp');
        
        % Inform the user
        fprintf(output,'%s\n',[pad('Number of slices in the interpolated control sequence',60) ...
            int2str(sum(ctrl_sys.cavity_n_interp))]);
        
    else
        if isfield(ctrl_param,'cavity_n_interp')
            ctrl_param=rmfield(ctrl_param,'cavity_n_interp');
        end
        if isfield(ctrl_param,'cavity_function')
            ctrl_param=rmfield(ctrl_param,'cavity_function');
        end
        
        ctrl_sys.cavity_n_interp=1+zeros(1,ctrl_param.pulse_nsteps);
    end
    
    
elseif isfield(ctrl_param,'cavity_n_interp')
    
    ctrl_param=rmfield(ctrl_param,'cavity_n_interp');
    
    % Default cavity decay rate is zero
    ctrl_sys.cavity_decay_rate=0;
    ctrl_sys.cavity_function=[];
    ctrl_sys.cavity_n_interp=1+zeros(1,ctrl_sys.pulse_nsteps);
    
else
    
    % Default cavity decay rate is zero
    ctrl_sys.cavity_decay_rate=0;
    ctrl_sys.cavity_function=[];
    ctrl_sys.cavity_n_interp=1+zeros(1,ctrl_sys.pulse_nsteps);
    
end

% use split operator method for interaction hamiltonian
if isfield(ctrl_param,'prop_ind') && size(ctrl_param.prop_ind,2)==2
    
    % Absorb the specification
    ctrl_sys.prop_ind=ctrl_param.prop_ind;
    ctrl_param=rmfield(ctrl_param,'prop_ind');
    
    % inform the user
    fprintf(output,'%s\n',[pad('Sparse propagator index store',60),'true']);

elseif contains(func2str(ctrl_sys.optimcon_fun),{'escalade','qoala'}) && isfield(ctrl_param,'approximation') && strcmp(ctrl_param.approximation,'no_coupling')

    if strcmp(ctrl_sys.space,'liouville') && strcmp(ctrl_sys.basis,'sphten')
        % NOTE: remove the unit state because it is not needed
        pform_rows=kron(3*(1:ctrl_sys.nspins)-2,ones(1,9))+kron(ones(1,3*ctrl_sys.nspins),[0 1 2]); pform_rows=pform_rows(:);
        pform_cols=kron((1:3*ctrl_sys.nspins),[1 1 1]); pform_cols=pform_cols(:);
        ctrl_sys.prop_ind=[pform_rows pform_cols];
        ctrl_sys.prop_grad_ind={ctrl_sys.prop_ind};
    elseif strcmp(ctrl_sys.space,'liouville') && strcmp(ctrl_sys.basis,'zeeman')
        error('propagator indices not yet coded for the zeeman basis in a liouville space')
    elseif strcmp(ctrl_sys.space,'hilbert') && strcmp(ctrl_sys.basis,'zeeman')
        pform_rows=kron(2*(1:ctrl_sys.nspins)-1,ones(1,4))+kron(ones(1,2*ctrl_sys.nspins),[0 1]); pform_rows=pform_rows(:);
        pform_cols=kron((1:2*ctrl_sys.nspins),[1 1]); pform_cols=pform_cols(:);
        ctrl_sys.prop_ind=[pform_rows pform_cols];
        ctrl_sys.prop_grad_ind={ctrl_sys.prop_ind};
    else
        error(['propagator indices not yet coded for the ' ctrl_sys.basis ' basis in a ' ctrl_sys.space ' space'])
    end

    % inform the user
    fprintf(output,'%s\n',[pad('Sparse propagator index store',60),'true']);
    
elseif contains(func2str(ctrl_sys.optimcon_fun),{'escalade','qoala'})

    ctrl_sys.prop_ind=cell2mat(propagator_ind(ctrl_sys.space,ctrl_sys.basis,1));
    ctrl_sys.prop_grad_ind=propagator_ind(ctrl_sys.space,ctrl_sys.basis,ctrl_sys.nspins);
    
    % inform the user
    fprintf(output,'%s\n',[pad('Sparse propagator index store',60),'true']);
    
else
    ctrl_sys.prop_ind=[];
    ctrl_sys.prop_grad_ind={ctrl_sys.prop_ind};
    fprintf(output,'%s\n',[pad('Sparse propagator index store',60),'false']);
    
end


% which controls (cols) affect which spins (rows)
if isfield(ctrl_param,'spin_control')
    
    % Absorb the specification
    ctrl_sys.spin_control=logical(ctrl_param.spin_control);
    ctrl_param=rmfield(ctrl_param,'spin_control');
    
elseif contains(func2str(ctrl_sys.optimcon_fun),{'escalade','qoala'})

    % not yet automated -- but could be
    error('need to provide an array instructing which operators (cols) affect which spins (rows)');
else
    ctrl_sys.spin_control=[];
end

if isfield(ctrl_param,'gradops')
    
    % Absorb the specification
    ctrl_sys.gradops=ctrl_param.gradops;
    ctrl_param=rmfield(ctrl_param,'gradops');
    
else
    ctrl_sys.gradops='full';
end


% use split operator method for interaction hamiltonian HERE TO DO
if isfield(ctrl_param,'pauli_operators')
    
    % Absorb the specification
    ctrl_sys.pauli_operators=ctrl_param.pauli_operators;
    ctrl_param=rmfield(ctrl_param,'pauli_operators');

elseif contains(func2str(ctrl_sys.optimcon_fun),{'escalade','qoala'})

    if strcmp(ctrl_sys.space,'liouville') && strcmp(ctrl_sys.basis,'sphten')

        L1x=[ 0 +1  0; +1  0 +1;  0 +1  0]; L1x=(1/sqrt(2))*L1x;
        L1y=[ 0 -1  0; +1  0 -1;  0 +1  0]; L1y=(1i/sqrt(2))*L1y;
        L1z=[+1  0  0;  0  0  0;  0  0 -1];

    elseif strcmp(ctrl_sys.space,'liouville') && strcmp(ctrl_sys.basis,'zeeman')

        L1x=[ 0 +1 -1 0; +1  0  0 -1; -1  0  0 +1;  0 -1 +1  0]; L1x=( 1/2)*L1x;
        L1y=[ 0 -1 -1 0; +1  0  0 -1; +1  0  0 -1;  0 +1 +1  0]; L1y=(1i/2)*L1y;
        L1z=[ 0  0  0 0;  0 -1  0  0;  0  0 +1  0;  0  0  0  0];

    elseif strcmp(ctrl_sys.space,'hilbert') && strcmp(ctrl_sys.basis,'zeeman')

        L1x=[ 0 +1; +1  0]; L1x=( 1/2)*L1x;
        L1y=[ 0 -1; +1  0]; L1y=(1i/2)*L1y;
        L1z=[+1  0;  0 -1]; L1z=( 1/2)*L1z;

    else
        error(['single spin operators not coded for the ' ctrl_sys.basis ' basis in a ' ctrl_sys.space ' space formalism.'])

    end

    % only coded for spin-1/2 at the moment
    for n=1:length(unique_mults)
        if unique_mults(n)~=2, error('only coded for spin 1/2'); end
    end

    if strcmp(ctrl_sys.gradops,'sparse')

        % pre-allocate the Pauli operator cell
        ctrl_sys.pauli_operators=cell(max(unique_mults)-1,3);
        
        for n=1:length(unique_mults)

            if unique_mults(n)~=2, error('only coded for spin 1/2'); end
            
            ctrl_sys.pauli_operators(unique_mults(n)-1,:)={L1x,L1y,L1z};
        end  
        
    elseif strcmp(ctrl_sys.gradops,'full') || strcmp(ctrl_sys.gradops,'sparse')
        
        % pre-allocate the Pauli operator cell
        ctrl_sys.pauli_operators=cell(size(ctrl_sys.spin_control,1),3*size(ctrl_sys.spin_control,2)/2);
        zmat=zeros(ctrl_sys.nspins,ctrl_sys.nspins);

        if strcmp(ctrl_sys.basis,'sphten')
            warning('[THIS MIGHT NOT BE CODED PROPERLY -- includind the unit state]')
                for k=1:size(ctrl_sys.spin_control,2)/2
                    for n=1:ctrl_sys.nspins
                        if ctrl_sys.spin_control(n,2*k)
                            if strcmp(ctrl_sys.gradops,'full')
                                mask_mat=zmat; mask_mat(ctrl_sys.nspins-n+1,ctrl_sys.nspins-n+1)=1;
                                ctrl_sys.pauli_operators{n,3*k-2}=sparse(blkdiag(0,kron(mask_mat,L1x)));
                                ctrl_sys.pauli_operators{n,3*k-1}=sparse(blkdiag(0,kron(mask_mat,L1y)));
                                ctrl_sys.pauli_operators{n,3*k}=sparse(blkdiag(0,kron(mask_mat,L1z)));
                            elseif strcmp(ctrl_sys.gradops,'sparse')
                                ctrl_sys.pauli_operators{n,3*k-2}=L1x;
                                ctrl_sys.pauli_operators{n,3*k-1}=L1y;
                                ctrl_sys.pauli_operators{n,3*k}=L1z;
                            end
                        end

                    end
                end

        elseif strcmp(ctrl_sys.basis,'zeeman')

                for k=1:size(ctrl_sys.spin_control,2)/2
                    for n=1:ctrl_sys.nspins
                        if ctrl_sys.spin_control(n,2*k)
                            if strcmp(ctrl_sys.gradops,'full')
                                mask_mat=zmat; mask_mat(ctrl_sys.nspins-n+1,ctrl_sys.nspins-n+1)=1;
                                ctrl_sys.pauli_operators{n,3*k-2}=sparse(kron(mask_mat,L1x));
                                ctrl_sys.pauli_operators{n,3*k-1}=sparse(kron(mask_mat,L1y));
                                ctrl_sys.pauli_operators{n,3*k}=sparse(kron(mask_mat,L1z));
                            elseif strcmp(ctrl_sys.gradops,'sparse')
                                ctrl_sys.pauli_operators{n,3*k-2}=L1x;
                                ctrl_sys.pauli_operators{n,3*k-1}=L1y;
                                ctrl_sys.pauli_operators{n,3*k}=L1z;
                            end
                        end

                    end
                end

        else
                error('only coded for zeeman or spherical tensor basies')

        end
    end
else
    ctrl_sys.pauli_operators={};

end

% Process offsets
if isfield(ctrl_param,'offsets') && isfield(ctrl_param,'offset_operator') &&...
        numel(ctrl_param.offsets)==numel(ctrl_param.offset_operator)
    
    % Absorb offsets and operators
    ctrl_sys.offsets=ctrl_param.offsets; 
    ctrl_param=rmfield(ctrl_param,'offsets');
    ctrl_sys.offset_operator=ctrl_param.offset_operator;
    ctrl_param=rmfield(ctrl_param,'offset_operator');
    
    % Inform the user and clean up the operators
    for n=1:numel(ctrl_sys.offset_operator)
        
        fprintf(output,'%s\n',[pad(['Number of systems in offset ensemble ' int2str(n)],60) ...
                            int2str(numel(ctrl_sys.offsets{n}))]);
        fprintf(output,'%s\n',[pad(['Max offset in the ensemble ' int2str(n) ', Hz'],60) ...
                            num2str(max(ctrl_sys.offsets{n}),'%+.9g')]);
        fprintf(output,'%s\n',[pad(['Min offset in the ensemble ' int2str(n) ', Hz'],60) ...
                            num2str(min(ctrl_sys.offsets{n}),'%+.9g')]);
    end
    
elseif isfield(ctrl_param,'offset_operator')

    % Default is no offsets
    ctrl_sys.offsets={0};
    ctrl_sys.offset_operator=ctrl_param.offset_operator;
    ctrl_param=rmfield(ctrl_param,'offset_operator');

else
    
    % Default is no offsets
    ctrl_sys.offsets={0};
    ctrl_sys.offset_operator={sparse(0)};
    
end

% symmetric pulses
if isfield(ctrl_param,'pulse_syms') && numel(ctrl_param.pulse_syms)==2 && sum(abs(ctrl_param.pulse_syms))==2
    
    % Absorb the specification
    ctrl_sys.pulse_syms=ctrl_param.pulse_syms;
    ctrl_param=rmfield(ctrl_param,'pulse_syms');
    
    % inform the user
    fprintf(output,'%s\n',[pad('Symmetric pulse constraint',60) 'true']);
    
else
    if isfield(ctrl_param,'pulse_syms')
        ctrl_param=rmfield(ctrl_param,'pulse_syms');
    end
    ctrl_sys.pulse_syms=[0 0];
end

if isfield(ctrl_param,'step_method')
    
    % Absorb expm step method
    ctrl_sys.step_method=ctrl_param.step_method;
    ctrl_param=rmfield(ctrl_param,'step_method');
    
else
    ctrl_sys.step_method='krylov';
end

if isfield(ctrl_param,'auxmat_method')
    
    % Absorb expm auxiliary matrix method
    ctrl_sys.auxmat_method=ctrl_param.auxmat_method;
    ctrl_param=rmfield(ctrl_param,'auxmat_method');
    
else
    ctrl_sys.auxmat_method='taylor';
end

% Optimisation method
if ~isfield(ctrl_param,'method')        
    
    % Default is LBFGS
    ctrl_sys.method='lbfgs';
    
else    
    
    % Absorb the method
    ctrl_sys.method=ctrl_param.method;
    ctrl_param=rmfield(ctrl_param,'method');
    
end

% Inform the user
fprintf(output,'%s\n',[pad('Optimisation method',60) pad(ctrl_sys.method,20)]);

% Maximum optimisation iterations
if isfield(ctrl_param,'max_iter')
    
    % Absorb the specification
    ctrl_sys.max_iter=ctrl_param.max_iter;
    ctrl_param=rmfield(ctrl_param,'max_iter');
    
else
    
    % Default is 100
    ctrl_sys.max_iter=100;
    
end

% Inform the user
fprintf(output,'%s\n',[pad('Maximum number of iterations',60) ...
                    pad(num2str(ctrl_sys.max_iter),20)]);

% Termination tolerance on the fidelity
if isfield(ctrl_param,'tol_f')
    
    % Absorb the spec
    ctrl_sys.tol_f=ctrl_param.tol_f;
    ctrl_param=rmfield(ctrl_param,'tol_f');
    
else
    
    % Default is 0.1%
    ctrl_sys.tol_f=inf;

end

% Inform the user
fprintf(output,'%s\n',[pad('Termination tolerance on fidelity',60) ...
                    pad(num2str(ctrl_sys.tol_f,'%0.8g'),20)]);

% Termination tolerance on the step norm
if isfield(ctrl_param,'tol_x')
    
    % Absorb the spec
    ctrl_sys.tol_x=ctrl_param.tol_x;
    ctrl_param=rmfield(ctrl_param,'tol_x');
    
else
    
    % Default is 0.1%
    ctrl_sys.tol_x=1e-3;

end

% Inform the user
fprintf(output,'%s\n',[pad('Termination tolerance on |delta_x|',60) ...
                    pad(num2str(ctrl_sys.tol_x,'%0.8g'),20)]);

% Termination tolerance on the gradient norm
if isfield(ctrl_param,'tol_g')
    
    % Absorb the spec
    ctrl_sys.tol_g=ctrl_param.tol_g;
    ctrl_param=rmfield(ctrl_param,'tol_g');
    
else
    
    % Default is 1e-6
    ctrl_sys.tol_g=1e-6;
    
end

% Inform the user
fprintf(output,'%s\n',[pad('Termination tolerance on |grad(x)|',60) ...
                    pad(num2str(ctrl_sys.tol_g,'%0.8g'),20)]);

% Set up LBFGS history
if strcmp(ctrl_sys.method,'lbfgs')
    
    % Decide history length
    if isfield(ctrl_param,'n_grads')
        
        % Absorb the spec
        ctrl_sys.n_grads=ctrl_param.n_grads;
        ctrl_param=rmfield(ctrl_param,'n_grads');
        
    else
        
        % Default is 20 gradients
        ctrl_sys.n_grads=20;
        
    end
    
    % Inform the user
    fprintf(output,'%s\n',[pad('Number of gradients in LBFGS history',60) ...
                        pad(num2str(ctrl_sys.n_grads),20)]);

elseif isfield(ctrl_param,'n_grads')
    
    % remove because it's not required
    ctrl_param=rmfield(ctrl_param,'n_grads');   

    % Inform the user
    fprintf(output,'%s\n',[pad(['WARNING: LBFGS history not used in ' ctrl_sys.method],60) ...
                        pad('(removed)',20)]);
end

% Linsearch parameter
if isfield(ctrl_param,'ls_c1')
    
    % Absorb the spec
    ctrl_sys.ls_c1=ctrl_param.ls_c1;
    ctrl_param=rmfield(ctrl_param,'ls_c1');

    % inform the user
    fprintf(output,'%s\n',[pad('Linesearch: Sufficient decrease condition, c1',60) ...
                    pad(num2str(ctrl_sys.ls_c1,'%0.8g'),20)]);
    
else
    
    % Default is 1e-6
    ctrl_sys.ls_c1=1e-2;
    
end

% Linsearch parameter
if isfield(ctrl_param,'ls_c2')
    
    % Absorb the spec
    ctrl_sys.ls_c2=ctrl_param.ls_c2;
    ctrl_param=rmfield(ctrl_param,'ls_c2');

    % inform the user
    fprintf(output,'%s\n',[pad('Linesearch: Curvature condition on grad, c2',60) ...
                    pad(num2str(ctrl_sys.ls_c2,'%0.8g'),20)]);
    
else
    
    % Default is 1e-6
    ctrl_sys.ls_c2=0.9;
    
end

% Linsearch parameter
if isfield(ctrl_param,'ls_tau1')
    
    % Absorb the spec
    ctrl_sys.ls_tau1=ctrl_param.ls_tau1;
    ctrl_param=rmfield(ctrl_param,'ls_tau1');

    % inform the user
    fprintf(output,'%s\n',[pad('Linesearch: Bracket expansion factor, tau1',60) ...
                    pad(num2str(ctrl_sys.ls_tau1,'%0.8g'),20)]);
    
else
    
    % Default is 1e-6
    ctrl_sys.ls_tau1=3;
    
end

% Linsearch parameter
if isfield(ctrl_param,'ls_tau2')
    
    % Absorb the spec
    ctrl_sys.ls_tau2=ctrl_param.ls_tau2;
    ctrl_param=rmfield(ctrl_param,'ls_tau2');

    % inform the user
    fprintf(output,'%s\n',[pad('Linesearch: Left section contraction, tau2',60) ...
                    pad(num2str(ctrl_sys.ls_tau2,'%0.8g'),20)]);
    
else
    
    % Default is 1e-6
    ctrl_sys.ls_tau2=0.1;
    
end

% Linsearch parameter
if isfield(ctrl_param,'ls_tau3')
    
    % Absorb the spec
    ctrl_sys.ls_tau3=ctrl_param.ls_tau3;
    ctrl_param=rmfield(ctrl_param,'ls_tau3');

    % inform the user
    fprintf(output,'%s\n',[pad('Linesearch: Right section contraction, tau3',60) ...
        pad(num2str(ctrl_sys.ls_tau3,'%0.8g'),20)]);

else

    % Default is 1e-6
    ctrl_sys.ls_tau3=0.5;

end

%if strcmp(ctrl_sys.method,'newton-raphson')
    % Max regularisation iterations - RFO
    if isfield(ctrl_param,'reg_max_iter')

        % Absorb the spec
        ctrl_sys.reg_max_iter=ctrl_param.reg_max_iter;
        ctrl_param=rmfield(ctrl_param,'reg_max_iter');

    else

        % Default is 1e-6
        ctrl_sys.reg_max_iter=2500;

    end

    % inform the user
    fprintf(output,'%s\n',[pad('Maximum Hessian regularisation iterations - RFO',60) ...
        pad(int2str(ctrl_sys.reg_max_iter),20)]);

    % RFO scaling factor - RFO
    if isfield(ctrl_param,'reg_alpha')

        % Absorb the spec
        ctrl_sys.reg_alpha=ctrl_param.reg_alpha;
        ctrl_param=rmfield(ctrl_param,'reg_alpha');

    else

        % Default is 1e-6
        ctrl_sys.reg_alpha=1;

    end

    % inform the user
    fprintf(output,'%s\n',[pad('Regularisation scaling factor - RFO',60) ...
        pad(num2str(ctrl_sys.reg_alpha,'%0.8g'),20)]);

    % RFO conditioning multiplier - RFO
    if isfield(ctrl_param,'reg_phi')

        % Absorb the spec
        ctrl_sys.reg_phi=ctrl_param.reg_phi;
        ctrl_param=rmfield(ctrl_param,'reg_phi');

    else

        % Default is 1e-6
        ctrl_sys.reg_phi=0.9;

    end

    % inform the user
    fprintf(output,'%s\n',[pad('Regularisation conditioning multiplier - RFO',60) ...
        pad(num2str(ctrl_sys.reg_phi,'%0.8g'),20)]);

    % Max. condition number - RFO
    if isfield(ctrl_param,'reg_max_cond')

        % Absorb the spec
        ctrl_sys.reg_max_cond=ctrl_param.reg_max_cond;
        ctrl_param=rmfield(ctrl_param,'reg_max_cond');

    else

        % Default is 1e-6
        ctrl_sys.reg_max_cond=1e4;

    end

    % inform the user
    fprintf(output,'%s\n',[pad('Maximum condition number - RFO',60) ...
        pad(num2str(ctrl_sys.reg_max_cond,'%0.8g'),20)]);
%end

% Accept pulse sequence parameters without question
if isfield(ctrl_param,'parameters')
    ctrl_sys.parameters=ctrl_param.parameters;
    ctrl_param=rmfield(ctrl_param,'parameters');
    fprintf(output,'%s\n','user parameters received without parsing.');
end

mem_req=0;
fprintf(output,'%s\n',pad('Calculating approximate memory footprint...',60));
% estimate the mem requirement for a (full) trajectory
mem_traj=((1+sum(ctrl_sys.cavity_n_interp))*8+16*ctrl_sys.dim*(1+sum(ctrl_sys.cavity_n_interp)));
[val,p]=findprefix(mem_traj,1024);
fprintf(output,'%s\n',[pad('... trajectory array',60) [pad(num2str(val,'%7.3f'),7,'right') ' ' p 'B']]);

mem_req=mem_req+mem_traj;

if contains(func2str(ctrl_sys.optimcon_fun),{'escalade','qoala'})

    % number of spins
    if isfield(ctrl_sys.drift_sys,'spn_mults')
        M1=length(ctrl_sys.drift_sys.spn_mults);
    else
        M1=1;
    end
    if isfield(ctrl_sys.drift_sys,'ind_spn')
        M2=length(ctrl_sys.drift_sys.ind_spn);
    else
        M2=1;
    end
    try
        spin_sys.bas.formalism=spin_system.bas.formalism;
        spin_sys.control.spin_control=spin_control;
        spin_sys.control.gradops=ctrl_sys.gradops;
        [P,dP]=rodrigues(spin_sys,ctrl_sys.pauli_operators,randn(3,M1*sum(ctrl_sys.cavity_n_interp)*ctrl_sys.nspins),randn(1,M1*sum(ctrl_sys.cavity_n_interp)*ctrl_sys.nspins)); %#ok<ASGLU> 
        P=whos('P'); dP=whos('dP');

        [val,p]=findprefix(P.bytes,1024);
        fprintf(output,'%s\n',[pad('... propagator elements array',60) [pad(num2str(val,'%7.3f'),7,'right') ' ' p 'B']]);
        
        [val,p]=findprefix(dP.bytes,1024);
        fprintf(output,'%s\n',[pad('... propagator derivative elements array',60) [pad(num2str(val,'%7.3f'),7,'right') ' ' p 'B']]);

        mem_req=mem_req+P.bytes+dP.bytes;
    catch
        mem_Pelements=(9*16*M1*ctrl_sys.nspins*sum(ctrl_sys.cavity_n_interp));
        
        [val,p]=findprefix(mem_Pelements,1024);
        fprintf(output,'%s\n',[pad('... propagator elements array',60) [pad(num2str(val,'%7.3f'),7,'right') ' ' p 'B']]);

        mem_req=mem_req+mem_Pelements;
    end

    mem_split=((M2*sum(ctrl_sys.cavity_n_interp))*8+M2*16*ctrl_sys.dim*(sum(ctrl_sys.cavity_n_interp)));
        
    [val,p]=findprefix(mem_split,1024);
    fprintf(output,'%s\n',[pad('... split trajectory array',60) [pad(num2str(val,'%7.3f'),7,'right') ' ' p 'B']]);

    mem_req=mem_req+mem_split;

    if (10^ctrl_sys.nspins/(4^ctrl_sys.nspins)^2)>0.15%spin_system.tols.dense_matrix
        mem_prop=M1*16*(ctrl_sys.dim^2)*(sum(ctrl_sys.cavity_n_interp));
        
        [val,p]=findprefix(mem_prop,1024);
        fprintf(output,'%s\n',[pad(['... ' num2str(sum(ctrl_sys.cavity_n_interp)) ' propagators (full matrices)'],60) [pad(num2str(val,'%7.3f'),7,'right') ' ' p 'B']]);

    else

        mem_prop=M1*(8*(1+ctrl_sys.dim) + (8+16)*(10*ctrl_sys.nspins))*(sum(ctrl_sys.cavity_n_interp));

        [val,p]=findprefix(mem_prop,1024);
        fprintf(output,'%s\n',[pad(['... ' num2str(sum(ctrl_sys.cavity_n_interp)) ' propagators (sparse matrices, no coupling)'],60) [pad(num2str(val,'%7.3f'),7,'right') ' ' p 'B']]);

    end
    mem_req=mem_req+mem_prop;

    [val,p]=findprefix(mem_req,1024);
    fprintf(output,'%s\n',[pad(['=== Total memory footprint of @' func2str(ctrl_sys.optimcon_fun)],60) [pad(num2str(val,'%7.3f'),7,'right') ' ' p 'B']]);

end

% Catch unparsed fields
unparsed=fieldnames(ctrl_param);
if ~isempty(unparsed)
    for n=1:numel(unparsed)
        fprintf(output,'%s\n',[pad(['Unrecognised option - ' unparsed{n}],60) 'WARNING']);
    end
end

end

% Consistency enforcement
function grumble(ctrl_param)
if ~isfield(ctrl_param,'optimcon_fun') || ~isa(ctrl_param.optimcon_fun,'function_handle')
    error('must provide the function handle to an optimal control function in ctrl_param.optimcon_fun.')
end
if isfield(ctrl_param,'space') && ~ismember(ctrl_param.space,{'liouville','hilbert'})
    error('only coded for Hilbert or Liouville spaces');
end
if isfield(ctrl_param,'basis') && ~ismember(ctrl_param.basis,{'sphten','zeeman'})
    error('only coded for Zeeman or spherical tensor basis sets');
end
if ~isfield(ctrl_param,'nspins') && ~isfield(ctrl_param,'mults')
    error('need to provide either ctrl_param.nspins or ctrl_param.mults');
end
if isfield(ctrl_param,'fidelity')
    if ~ischar(ctrl_param.fidelity)
        error('control.fidelity must be a character string.');
    end
    if ~ismember(ctrl_param.fidelity,{'real','square'})
        error('control.fidelity can be ''real'' or ''square''.');
    end
end
if ~isfield(ctrl_param,'operators')
    error('control operators must be supplied in control.operators field.');
elseif isfield(ctrl_param,'optimcon_fun') && contains(func2str(ctrl_param.optimcon_fun),'escalade')
    if ~isfield(ctrl_param,'ctrl_axes')
        error('control axes must be supplied in control.axes field if an escalade method is used.');
    end
end
if (~iscell(ctrl_param.operators))||...
   (~all(cellfun(@ismatrix,ctrl_param.operators(:))))
    error('control.operators must be a cell array of matrices.');
end
if isfield(ctrl_param,'pulse_dt')
    if isfield(ctrl_param,'pulse_dur')
        error('control.pulse_dur cannot be specified simultaneously with control.pulse_dt');
    end
    if isfield(ctrl_param,'pulse_nsteps')
        error('control.pulse_nsteps cannot be specified simultaneously with control.pulse_dt');
    end
    if (~isnumeric(ctrl_param.pulse_dt))||(~isreal(ctrl_param.pulse_dt))||...
       (~iscolumn(ctrl_param.pulse_dt))||any(ctrl_param.pulse_dt(:)<=0)
        error('control.pulse_dt must be a column vector of positive real numbers.');
    end
else
    if ~isfield(ctrl_param,'pulse_dur')
        error('pulse duration must be specified in control.pulse_dur or control.pulse_dt field.');
    end
    if (~isnumeric(ctrl_param.pulse_dur))||(~isreal(ctrl_param.pulse_dur))||...
       (~isscalar(ctrl_param.pulse_dur))||(ctrl_param.pulse_dur<=0)
        error('control.pulse_dur must be a positive real scalar.');
    end
    if ~isfield(ctrl_param,'pulse_nsteps')
        error('number of pulse slices must be specified in control.pulse_nsteps or control.pulse_dt field.');
    end
    if (~isnumeric(ctrl_param.pulse_nsteps))||(~isreal(ctrl_param.pulse_nsteps))||...
       (~isscalar(ctrl_param.pulse_nsteps))||(ctrl_param.pulse_nsteps<1)
        error('control.pulse_nsteps must be a positive real integer.');
    end
end
if ~isfield(ctrl_param,'pulse_nsteps')
    nsteps=length(ctrl_param.pulse_dt);
else
    nsteps=ctrl_param.pulse_nsteps;
end
if isfield(ctrl_param,'drift_sys')
    if iscell(ctrl_param.drift_sys)
        for n=1:numel(ctrl_param.drift_sys)
            drift_check=ctrl_param.drift_sys{n};
            if ~isfield(drift_check,'drift') && (~isfield(drift_check,'interaction') || ~isfield(drift_check,'singlespin')) && ~isfield(drift_check,'offset')
                error(['control.drift_sys{' num2str(n) '} must contain at least either the field drift, interaction and/or singlespin, or offset'])
            end
            if isfield(ctrl_param,'optimcon_fun') && ~contains(func2str(ctrl_param.optimcon_fun),'escalade')
                if ~isfield(drift_check,'offset')
                    error(['control.drift_sys{' num2str(n) '}.offset field must be given for an escalade method'])
                end
                if iscell(drift_check.offset)
                    if numel(drift_check.offset)~=nsteps
                        error(['control.drift_sys{' num2str(n) '}.offset as a cell array must have elements equal to the number of time steps'])
                    end
                    for m=1:numel(drift_check.offset)
                        if ~isnumeric(drift_check.offset{m}) || ~isreal(drift_check.offset{m}) || ~iscolumn(drift_check.offset{m})
                            error(['control.drift_sys{' num2str(n) '}.offset{' num2str(m) '} must be a column vector of real numbers'])
                        end
                    end
                elseif ~isnumeric(drift_check.offset) || ~isreal(drift_check.offset) || ~iscolumn(drift_check.offset)
                    error(['control.drift_sys{' num2str(n) '}.offset must be a cell array (time dependence) or column vector of real numbers'])
                end
                if isfield(drift_check,'interaction')
                    if iscell(drift_check.interaction)
                        if numel(drift_check.interaction)~=nsteps
                            error(['control.drift_sys{' num2str(n) '}.interaction as a cell array must have elements equal to the number of time steps'])
                        end
                        for m=1:numel(drift_check.interaction)
                            if ~isnumeric(drift_check.interaction{m}) || ~ismatrix(drift_check.interaction) || size(drift_check.interaction{m},1)~=size(drift_check.interaction{m},2)
                                error(['control.drift_sys{' num2str(n) '}.interaction{' num2str(m) '} must be a square matrix of complex numbers'])
                            end
                        end
                    elseif ~isnumeric(drift_check.interaction) || ~ismatrix(drift_check.interaction) || size(drift_check.interaction,1)~=size(drift_check.interaction,2)
                        error(['control.drift_sys{' num2str(n) '}.interaction must be a cell array (time dependence) or square matrix of complex numbers'])
                    end
                end
                % adaptivity parameters
                if isfield(drift_check,'adapt_method')
                    if ~ismember(drift_check.adapt_method,{'gain','exact'})
                        error(['allowable methods of control.drift_sys{' num2str(n) '}.adapt_method are ''gain'' and ''exact''']);
                    end
                end
                if isfield(drift_check,'adapt_scale')
                    if ~all(ismember(drift_check.adapt_scale,{'none','infidelity','gradient'}))
                        error(['allowable methods of control.drift_sys{' num2str(n) '}.adapt_scale are ''none'',''infidelity'', or ''gradient''']);
                    end
                end
            else
                if ~isfield(drift_check,'drift') && (~isfield(drift_check,'interaction') || ~isfield(drift_check,'singlespin'))
                    error(['control.drift_sys{' num2str(n) '} must contain either the field drift, or the fields interaction and/or singlespin'])
                end
                if isfield(drift_check,'drift')
                    if iscell(drift_check.drift)
                        if numel(drift_check.drift)~=nsteps
                            if numel(drift_check.drift)~=1
                                error(['control.drift_sys{' num2str(n) '}.drift as a cell array must have elements equal to the number of time steps'])
                            end
                        end
                        for m=1:numel(drift_check.drift)
                            if ~isnumeric(drift_check.drift{m}) || ~ismatrix(drift_check.drift) || size(drift_check.drift{m},1)~=size(drift_check.drift{m},2)
                                error(['control.drift_sys{' num2str(n) '}.drift{' num2str(m) '} must be a square matrix of complex numbers'])
                            end
                        end
                    elseif ~isnumeric(drift_check.drift) || ~ismatrix(drift_check.drift) || size(drift_check.drift,1)~=size(drift_check.drift,2)
                        error(['control.drift_sys{' num2str(n) '}.drift must be a cell array (time dependence) or square matrix of complex numbers'])
                    end
                end
                if isfield(drift_check,'singlespin')
                    if iscell(drift_check.singlespin)
                        if numel(drift_check.singlespin)~=nsteps
                            error(['control.drift_sys{' num2str(n) '}.singlespin as a cell array must have elements equal to the number of time steps'])
                        end
                        for m=1:numel(drift_check.singlespin)
                            if ~isnumeric(drift_check.singlespin{m}) || ~ismatrix(drift_check.singlespin) || size(drift_check.singlespin{m},1)~=size(drift_check.singlespin{m},2)
                                error(['control.drift_sys{' num2str(n) '}.singlespin{' num2str(m) '} must be a square matrix of complex numbers'])
                            end
                        end
                    elseif ~isnumeric(drift_check.singlespin) || ~ismatrix(drift_check.singlespin) || size(drift_check.singlespin,1)~=size(drift_check.singlespin,2)
                        error(['control.drift_sys{' num2str(n) '}.singlespin must be a cell array (time dependence) or square matrix of complex numbers'])
                    end
                end
                if isfield(drift_check,'interaction')
                    if iscell(drift_check.interaction)
                        if numel(drift_check.interaction)~=nsteps
                            error(['control.drift_sys{' num2str(n) '}.interaction as a cell array must have elements equal to the number of time steps'])
                        end
                        for m=1:numel(drift_check.interaction)
                            if ~isnumeric(drift_check.interaction{m}) || ~ismatrix(drift_check.interaction) || size(drift_check.interaction{m},1)~=size(drift_check.interaction{m},2)
                                error(['control.drift_sys{' num2str(n) '}.interaction{' num2str(m) '} must be a square matrix of complex numbers'])
                            end
                        end
                    elseif ~isnumeric(drift_check.interaction) || ~ismatrix(drift_check.interaction) || size(drift_check.interaction,1)~=size(drift_check.interaction,2)
                        error(['control.drift_sys{' num2str(n) '}.interaction must be a cell array (time dependence) or square matrix of complex numbers'])
                    end
                end
                % adaptivity parameters
                if isfield(drift_check,'adapt_method')
                    if ~ismember(drift_check.adapt_method,{'gain','exact'})
                        error(['allowable methods of control.drift_sys{' num2str(n) '}.adapt_method are ''gain'' and ''exact''']);
                    end
                end
                if isfield(drift_check,'adapt_scale')
                    if ~all(ismember(drift_check.adapt_scale,{'none','infidelity','gradient'}))
                        error(['allowable methods of control.drift_sys{' num2str(n) '}.adapt_scale are ''none'',''infidelity'', or ''gradient''']);
                    end
                end
            end
        end
    else
        drift_check=ctrl_param.drift_sys;
        if ~isfield(drift_check,'drift') && (~isfield(drift_check,'interaction') || ~isfield(drift_check,'singlespin')) && ~isfield(drift_check,'offset')
            error('control.drift_sys must contain at least either the field drift, interaction and/or singlespin, or offset')
        end
        if isfield(ctrl_param,'optimcon_fun') && contains(func2str(ctrl_param.optimcon_fun),'escalade')
            if isfield(drift_check,'offset')
                if iscell(drift_check.offset)
                    for m=1:numel(drift_check.offset)
                        if ~isnumeric(drift_check.offset{m}) || ~isreal(drift_check.offset{m}) || ~iscolumn(drift_check.offset{m})
                            error(['control.drift_sys.offset{' num2str(m) '} must be a column vector of real numbers'])
                        end
                    end
                elseif ~isnumeric(drift_check.offset) || ~isreal(drift_check.offset) || ~iscolumn(drift_check.offset)
                    error('control.drift_sys.offset must be a cell array (time dependence) or column vector of real numbers')
                end
            end
            if isfield(drift_check,'interaction')
                if iscell(drift_check.interaction)
                    for m=1:numel(drift_check.interaction)
                        if ~isnumeric(drift_check.interaction{m}) || ~ismatrix(drift_check.interaction) || size(drift_check.interaction{m},1)~=size(drift_check.interaction{m},2)
                            error(['control.drift_sys.interaction{' num2str(m) '} must be a square matrix of complex numbers'])
                        end
                    end
                elseif ~isnumeric(drift_check.interaction) || ~ismatrix(drift_check.interaction) || size(drift_check.interaction,1)~=size(drift_check.interaction,2)
                    error('control.drift_sys.interaction must be a cell array (time dependence) or square matrix of complex numbers')
                end
            end
            % adaptivity parameters
            if isfield(drift_check,'adapt_method')
                if ~ismember(drift_check.adapt_method,{'gain','exact'})
                    error('allowable methods of control.drift_sys.adapt_method are ''gain'' and ''exact''');
                end
            end
            if isfield(drift_check,'adapt_scale')
                if ~all(ismember(drift_check.adapt_scale,{'none','infidelity','gradient'}))
                    error('allowable methods of control.drift_sys.adapt_scale are ''none'',''infidelity'', or ''gradient''');
                end
            end
        else
            if ~isfield(drift_check,'drift') && ~isfield(drift_check,'interaction') && ~isfield(drift_check,'singlespin')
                error('control.drift_sys must contain at least the field drift, interaction, or singlespin')
            end
            if isfield(drift_check,'drift')
                if iscell(drift_check.drift)
                    for m=1:numel(drift_check.drift)
                        if ~isnumeric(drift_check.drift{m}) || ~ismatrix(drift_check.drift) || size(drift_check.drift{m},1)~=size(drift_check.drift{m},2)
                            error(['control.drift_sys.drift{' num2str(m) '} must be a square matrix of complex numbers'])
                        end
                    end
                elseif ~isnumeric(drift_check.drift) || ~ismatrix(drift_check.drift) || size(drift_check.drift,1)~=size(drift_check.drift,2)
                    error('control.drift_sys.drift must be a cell array (time dependence) or square matrix of complex numbers')
                end
            end
            if isfield(drift_check,'singlespin')
                if iscell(drift_check.singlespin)
                    for m=1:numel(drift_check.singlespin)
                        if ~isnumeric(drift_check.singlespin{m}) || ~ismatrix(drift_check.singlespin) || size(drift_check.singlespin{m},1)~=size(drift_check.singlespin{m},2)
                            error(['control.drift_sys.singlespin{' num2str(m) '} must be a square matrix of complex numbers'])
                        end
                    end
                elseif ~isnumeric(drift_check.singlespin) || ~ismatrix(drift_check.singlespin) || size(drift_check.singlespin,1)~=size(drift_check.singlespin,2)
                    error('control.drift_sys.singlespin must be a cell array (time dependence) or square matrix of complex numbers')
                end
            end
            if isfield(drift_check,'interaction')
                if iscell(drift_check.interaction)
                    for m=1:numel(drift_check.interaction)
                        if ~isnumeric(drift_check.interaction{m}) || ~ismatrix(drift_check.interaction) || size(drift_check.interaction{m},1)~=size(drift_check.interaction{m},2)
                            error(['control.drift_sys.interaction{' num2str(m) '} must be a square matrix of complex numbers'])
                        end
                    end
                elseif ~isnumeric(drift_check.interaction) || ~ismatrix(drift_check.interaction) || size(drift_check.interaction,1)~=size(drift_check.interaction,2)
                    error('control.drift_sys.interaction must be a cell array (time dependence) or square matrix of complex numbers')
                end
            end
            % adaptivity parameters
            if isfield(drift_check,'adapt_method')
                if ~ismember(drift_check.adapt_method,{'gain','exact'})
                    error('allowable methods of control.drift_sys.adapt_method are ''gain'' and ''exact''');
                end
            end
            if isfield(drift_check,'adapt_scale')
                if ~all(ismember(drift_check.adapt_scale,{'none','infidelity','gradient'}))
                    error('allowable methods of control.drift_sys.adapt_scale are ''none'',''infidelity'', or ''gradient''');
                end
            end
        end
    end
else
    error('information on the drift terms must be supplied in control.drift_sys field.');
end
if ~isfield(ctrl_param,'rho_targ') && ~isfield(ctrl_param,'prop_targ')
    error('targets must be supplied in control.rho_targ or control.prop_targ fields.');
end
if isfield(ctrl_param,'prop_targ') && ~iscell(ctrl_param.prop_targ)
    error('control.rho_targ must be a cell array of vectors.');
elseif isfield(ctrl_param,'prop_targ')
    for n=1:numel(ctrl_param.prop_targ)
        if (~isnumeric(ctrl_param.prop_targ{n}))||(~ismatrix(ctrl_param.prop_targ{n}))
            error('control.prop_targ must be a cell array of matrices.');
        end
    end
else
    if ~isfield(ctrl_param,'rho_init')
        error('initial states must be supplied in control.rho_init field.');
    end
    if ~iscell(ctrl_param.rho_init)
        error('control.rho_init must be a cell array of vectors.');
    end
    for n=1:numel(ctrl_param.rho_init)
        if strcmp(ctrl_param.space,'liouville')
            if (~isnumeric(ctrl_param.rho_init{n}))||(~iscolumn(ctrl_param.rho_init{n}))
                error('control.rho_init must be a cell array of column vectors.');
            end
        elseif strcmp(ctrl_param.space,'hilbert')
            if (~isnumeric(ctrl_param.rho_init{n}))||(~ismatrix(ctrl_param.rho_init{n}))||(size(ctrl_param.rho_init{n},1)~=size(ctrl_param.rho_init{n},2))
                error('control.rho_init must be a cell array of square matrices.');
            end
        end
    end
    if ~isfield(ctrl_param,'rho_targ')
        error('target states must be supplied in control.rho_targ field.');
    end
    if ~iscell(ctrl_param.rho_targ)
        error('control.rho_targ must be a cell array of vectors.');
    end
    for n=1:numel(ctrl_param.rho_targ)
        if strcmp(ctrl_param.space,'liouville')
            if (~isnumeric(ctrl_param.rho_targ{n}))||(~iscolumn(ctrl_param.rho_targ{n}))
                error('control.rho_targ must be a cell array of column vectors.');
            end
        elseif strcmp(ctrl_param.space,'hilbert')
            if (~isnumeric(ctrl_param.rho_targ{n}))||(~ismatrix(ctrl_param.rho_targ{n}))||(size(ctrl_param.rho_targ{n},1)~=size(ctrl_param.rho_targ{n},2))
                error('control.rho_targ must be a cell array of square matrices.');
            end
        end
    end
    if numel(ctrl_param.rho_targ)~=numel(ctrl_param.rho_init)
        error('control.rho_targ must have the same number of vectors as control.rho_init');
    end
end

if ~isfield(ctrl_param,'pwr_levels')
    error('power levels must be specified in control.pwr_levels field.');
end
if iscell(ctrl_param.pwr_levels) 
    if numel(ctrl_param.pwr_levels)~=numel(ctrl_param.operators)
        error('a cell array of power levels must have a number of elements equal to the number of control channels')
    end
    for n=1:numel(ctrl_param.pwr_levels)
        if (~isnumeric(ctrl_param.pwr_levels{n}))||(~isreal(ctrl_param.pwr_levels{n}))||...
        (~isrow(ctrl_param.pwr_levels{n}))||any(ctrl_param.pwr_levels{n}(:)<0)
            error(['power level array ' int2str(n) ' must be a row vector of positive real numbers'])
        end
    end
elseif (~isnumeric(ctrl_param.pwr_levels))||(~isreal(ctrl_param.pwr_levels))||...
   ((numel(ctrl_param.pwr_levels)~=1)&&(size(ctrl_param.pwr_levels,2)~=numel(ctrl_param.operators)))||any(ctrl_param.pwr_levels(:)<=0)
    error('if not a cell array, power levels must be a single number, or an array of K column vectors of positive real numbers, where K is the number of control channels.');
end
if isfield(ctrl_param,'penalties')
    if (~iscell(ctrl_param.penalties))||any(~cellfun(@ischar,ctrl_param.penalties(:)))
        error('control.penalties must be a cell array of character strings.');
    end
    if any(~ismember(ctrl_param.penalties(:),{'none','NS','SNS','DNS','SNSA','ADIAB'}))
        error('elements of control.penalties can be ''none'',''NS'', ''SNS'', ''SNSA'',''ADIAB'', or ''DNS''');
    end
end
if isfield(ctrl_param,'p_weights')
    if (~isnumeric(ctrl_param.p_weights))||(~isreal(ctrl_param.p_weights))||...
            (~isrow(ctrl_param.p_weights))||any(ctrl_param.p_weights(:)<=0)
        error('control.p_weights must be a row vector of positive real numbers.');
    end
    if ~isfield(ctrl_param,'penalties')
        error('penalties must be specified in control.penalties field.');
    end
end
if isfield(ctrl_param,'u_bound')
    if (~isnumeric(ctrl_param.u_bound))||(~isreal(ctrl_param.u_bound))||...
            (~isscalar(ctrl_param.u_bound) && numel(ctrl_param.u_bound)~=numel(ctrl_param.operators)*ctrl_param.pulse_nsteps)
        error('control.u_bound must be a real scalar or array.');
    end
    if ~isfield(ctrl_param,'l_bound')
        error('lower bound must be specified in control.l_bound field.');
    end
end
if isfield(ctrl_param,'l_bound')
    if (~isnumeric(ctrl_param.l_bound))||(~isreal(ctrl_param.l_bound))||...
            (~isscalar(ctrl_param.l_bound) && numel(ctrl_param.u_bound)~=numel(ctrl_param.operators)*ctrl_param.pulse_nsteps)
        error('control.l_bound must be a real scalar or array.');
    end
    if ~isfield(ctrl_param,'u_bound')
        error('upper bound must be specified in control.u_bound field.');
    end
end
if isfield(ctrl_param,'cavity_decay_rate')
    if (~isnumeric(ctrl_param.cavity_decay_rate)) ||...
            (~isreal(ctrl_param.cavity_decay_rate)) ||...
            ctrl_param.cavity_decay_rate<0
        error('control.cavity_decay_rate must be a real positive number');
    end
    if ctrl_param.cavity_decay_rate>0
        if ~isfield(ctrl_param,'cavity_n_interp')
            error('interpolation points control.cavity_n_interp must be supplied.')
        else
            if ~isreal(ctrl_param.cavity_n_interp) ||...
                    any(mod(ctrl_param.cavity_n_interp,1)~=0) ||...
                    any(ctrl_param.cavity_n_interp<=0)
                error('interpolation points control.cavity_n_interp must be an array of positive integer.');
            end
        end
        if ~isfield(ctrl_param,'cavity_function') || ~isa(ctrl_param.cavity_function,'function_handle')
            error('cavity_function must be supplied as a function handle.');
        end
    end
end
if isfield(ctrl_param,'sum_up')
    if ctrl_param.sum_up~=0 && ctrl_param.sum_up~=1
        error('control.sum_up must be either 0 or 1.');
    end
end
end

function str=helper(ctrl_param)

A = readmatrix('../readme.tex','FileType','text','Delimiter','\n','OutputType','string');
f = contains(A,'\item \mcode{param.');
A = A(f);

if ~isempty(ctrl_param)
    f = contains(A,['\item \mcode{param.' ctrl_param '}']);
    A=A(f);
end
% clean up
if ~isempty(A)
A = replace(A,'\mcode','');
A = replace(A,'\item ','');
A = replace(A,'\dots','...');
A = replace(A,'$','');
A = replace(A,'\leq','<=');
A = replace(A,'\geq','>=');
A = replace(A,'\neq','~=');
A = replace(A,'\textbf','');
A = replace(A,'\emph','');
A = replace(A,'\{','<lcurl>');
A = replace(A,'\}','<rcurl>');
A = replace(A,'{','');
A = replace(A,'}','');
A = replace(A,'<lcurl>','{');
A = replace(A,'<rcurl>','}');
A = replace(A,'\hfill','');
A = split(A,'\\ ');
A = replace(A,'\\','');

C=cell(1); C{1}=char(['<strong>' A(1) '</strong>']); C{1}=string(C{1});
for n=2:length(A)
    B=char(A(n));
    while ~isempty(B)
        f=isspace(B);
        if length(f)>60
            if sum(f(60:end))==0
                C{end+1}=['    ' B]; %#ok<AGROW> 
                B=[];
            else
                ind=59+find(f(60:end), 1, 'first');
                C{end+1}=['    ' B(1:ind-1)]; %#ok<AGROW> 
                B(1:ind)=[];
            end
        else
            C{end+1}=['    ' B]; %#ok<AGROW> 
            B=[];
        end
    end
    C{end+1}=''; %#ok<AGROW> 
end
else
    C{1}='';
    C{2}=char(['<strong>param.' ctrl_param '</strong>']); C{1}=string(C{1});
    C{3}='';
    C{4}=char('    Either this isn''t a parameter, or it isn''t in the readme.');
    C{5}='';
end
str=C;
end

function [out,prefix,prefixfull]=findprefix(in,Bas)
% NOTE: small problem with rounding and giving correct significant figures, probably needs attention at some point
if nargin<2
    Bas=1000;
end
if     round(in*(Bas^8))/(Bas^8)>=(Bas^8),          prefix='Y'; prefixfull='yotta'; out=in/(Bas^8);
elseif round(in*(Bas^7))/(Bas^7)>=(Bas^7),          prefix='Z'; prefixfull='zetta'; out=in/(Bas^7);
elseif round(in*(Bas^6))/(Bas^6)>=(Bas^6),          prefix='E'; prefixfull='exa';   out=in/(Bas^6);
elseif round(in*(Bas^5))/(Bas^5)>=(Bas^5),          prefix='P'; prefixfull='peta';  out=in/(Bas^5);
elseif round(in*(Bas^4))/(Bas^4)>=(Bas^4),          prefix='T'; prefixfull='tera';  out=in/(Bas^4);
elseif round(in*(Bas^3))/(Bas^3)>=(Bas^3),          prefix='G'; prefixfull='giga';  out=in/(Bas^3);
elseif round(in*(Bas^2))/(Bas^2)>=(Bas^2),          prefix='M'; prefixfull='mega';  out=in/(Bas^2);
elseif round(in*(Bas^1))/(Bas^1)>=(Bas^1),          prefix='k'; prefixfull='kilo';  out=in/(Bas^1);
elseif round(in*(Bas^0))/(Bas^0)>=(Bas^0) || in==0, prefix='';  prefixfull='';      out=in/(Bas^0);
elseif round(in/(Bas^-1))*(Bas^-1)>=(Bas^-1),       prefix='m'; prefixfull='milli'; out=in/(Bas^-1);
elseif round(in/(Bas^-2))*(Bas^-2)>=(Bas^-2),       prefix='u'; prefixfull='micro'; out=in/(Bas^-2);
elseif round(in/(Bas^-3))*(Bas^-3)>=(Bas^-3),       prefix='n'; prefixfull='nano';  out=in/(Bas^-3);
elseif round(in/(Bas^-4))*(Bas^-4)>=(Bas^-4),       prefix='p'; prefixfull='pico';  out=in/(Bas^-4);
elseif round(in/(Bas^-5))*(Bas^-5)>=(Bas^-5),       prefix='f'; prefixfull='femto'; out=in/(Bas^-5);
elseif round(in/(Bas^-6))*(Bas^-6)>=(Bas^-6),       prefix='a'; prefixfull='atto';  out=in/(Bas^-6);
elseif round(in/(Bas^-7))*(Bas^-7)>=(Bas^-7),       prefix='z'; prefixfull='zepto'; out=in/(Bas^-7);
else,                                               prefix='y'; prefixfull='yocto'; out=in/(Bas^-8);
end
end
