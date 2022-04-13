% Finds a local minimum of a function of several variables using Newton
% and quasi-Newton algorithms. Syntax:
%
%          [x,data]=fmaxnewton(param,cost_function,guess)
%
% Parameters:
%
%    param         - data object that has been
%                    through optimconset.m function
%
%    cost_function - a function handle that takes the input
%                    the size of guess
%
%    guess         - the initial point of the optimisation
%
% Outputs:
%
%    x               - the final point of the optimisation
%
%    data.count.iter - iteration counter
%
%    data.count.fx   - function evaluation counter
%
%    data.count.gfx  - gradient evaluation counter
%
%    data.count.hfx  - Hessian evaluation counter
%
%    data.count.rfo  - RFO iteration counter
%
%    data.x_shape    - output of size(guess)
%
%    data.*          - further fields may be set by the
%                      objective functon
%
% Author:       David L. Goodwin 
% Dev. period:  01/2020--04/2022
% Updated:      13/04/22
%
% Contact:
%   david.goodwin@chem.ox.ac.uk
%   david.goodwin@partner.kit.edu
%
% NOTE: This function is based on fminlbfgs.m code from D. Kroon,
%       University of Twente (Nov 2010), also implemented in Spinach
%       (Goodwin, Kuprov) and should be referenced accordingly.

function [x,data,exitflag]=fmaxnewton(param,cost_function,guess)

% Check consistency
grumble(param,cost_function,guess);

report_struct=[0 0 0];
if isfield(param,'fidelity_chk') &&...
        isa(param.fidelity_chk,'function_handle')
    report_struct(1)=1;
end
if numel(param.penalties)==1 &&...
        strcmp(param.penalties{1},'none')
else
    report_struct(2)=1;
end
if isfield(param.drift_sys{1},'adapt')
    report_struct(3)=1;
end

% Initialise counters
data.count.iter=0; data.count.fx=0;
data.count.gfx=0;  data.count.hfx=0;
data.count.rfo=0;

% set fidelity store
data.fx_store=NaN*ones(param.max_iter+1,2+length(param.penalties));
data.fx_chk_store=NaN*ones(param.max_iter+1,1+length(param.penalties));
data.timer=NaN*ones(param.max_iter+1,4);

% Stretch the guess supplied by user
data.x_shape=size(guess); x=guess(:);

% Print the header
header(param,report_struct);

% start the timer
opt_time=tic;

if param.max_iter==0
    
    % initialise this timer
    data.timer(1,1)=0;
    
    % Get objective and gradient
    [data,param,fx]=objective(x,cost_function,data,param);
    
    % store the fidelities
    data.fx_store(1,:)=[0 data.fx_sep_pen];
    data.timer(1,2)=toc(opt_time);
    data.fx_chk_store(1,:)=data.fx_chk;
    
    % Report diagnostics to user
    itrep(param,fx,[],[],data,report_struct)
    
    n=0;
else
    
    % Start the iteration loop
    for n=1:param.max_iter
        
        % initialise this timer
        if n==1,    data.timer(1,1)=0;
        else,       data.timer(n+1,1)=n;
        end
        
        
        % Get the search direction
        switch param.method
            
            case 'lbfgs'
                
                if n==1
                    
                    % Get objective and gradient
                    [data,param,fx,g]=objective(x,cost_function,data,param);
                    
                    % store the fidelities
                    data.fx_store(1,:)=[0 data.fx_sep_pen];
                    data.fx_chk_store(1,:)=data.fx_chk;
                    data.timer(1,2)=toc(opt_time);
                    
                    % Start history arrays
                    old_x=x; dx_hist=[]; old_g=g; dg_hist=[];
                    
                    % Direction vector
                    dir=g;
                    
                    % Report diagnostics to user
                    itrep(param,fx,g,[],data,report_struct)
                    
                    % set the next timer
                    data.timer(n+1,1)=n;
                    
                else
                    
                    % Update the history of dx and dg
                    dx_hist=[x-old_x dx_hist]; old_x=x; %#ok<AGROW>
                    dg_hist=[g-old_g dg_hist]; old_g=g; %#ok<AGROW>
                    
                    % Truncate the history
                    if size(dx_hist,2)>param.n_grads
                        dx_hist=dx_hist(:,1:param.n_grads);
                    end
                    if size(dg_hist,2)>param.n_grads
                        dg_hist=dg_hist(:,1:param.n_grads);
                    end
                    
                    % Get the direction
                    dir=search_direction('lbfgs',dx_hist,dg_hist,g);
                    
                end
                
            case 'newton-raphson'
                
                % Get objective, gradient, and Hessian
                [data,param,fx,g,H]=objective(x,cost_function,data,param);
                
                if n==1
                    
                    % store the fidelities
                    data.fx_store(1,:)=[0 data.fx_sep_pen];
                    data.timer(1,2)=toc(opt_time);
                    data.fx_chk_store(1,:)=data.fx_chk;
                    
                    % Report diagnostics to user
                    itrep(param,fx,g,[],data,report_struct)
                    
                    % set the next timer
                    data.timer(n+1,1)=n;
                    
                end
                
                % ensure real gradient and Hessian, and symmetric Hessian
                H=real(H+H')/2; g=real(g);
                
                % Regularise the Hessian
                [H,data]=hessianreg(param,-H,g,data);
                
                % Get the search direction
                dir=H\g;
                
        end
        
        % Update iteration counter
        data.count.iter=data.count.iter+1;

        % Run line search, refresh objective and gradient
        [alpha,fx,g,exitflag,data,param]=fmaxlinesearch(param,cost_function,dir,x,fx,g,data);
        
        % store after linesearch
        data.fx_store(n+1,:)=[n data.fx_sep_pen];
        data.fx_chk_store(n+1,:)=data.fx_chk;
        data.timer(n+1,2)=toc(opt_time);
        
        % Report diagnostics to user
        itrep(param,fx,g,alpha,data,report_struct)
        
        % Update x
        x=x+alpha*dir;
        
        % Check termination conditions
        if norm(alpha*dir,1)<param.tol_x, exitflag=2;
        elseif isempty(g), exitflag=1;
        elseif norm(g,2)<param.tol_g, exitflag=1;
        elseif fx>param.tol_f, exitflag=3;
        end
        
        % Exit if necessary
        if exitflag, break; end
        
    end
    
end

% See if iteration count was exceeded
if isempty(n)||(n==param.max_iter), exitflag=0; end

% Fold back the waveform
x=reshape(x,data.x_shape);

% Print the footer
footer(param,exitflag,data,report_struct);

end

% Header printing function
function header(param,report_struct)

if report_struct(3)==0
    if sum(report_struct)==2
        report_optim(param,'======================================================================================================');
        report_optim(param,'Iter  #f   #g   #H   #R    fidelity      check         penalise      total       alpha      |grad|    ');
        report_optim(param,'------------------------------------------------------------------------------------------------------');
    elseif sum(report_struct)==0
        report_optim(param,'==============================================================');
        report_optim(param,'Iter  #f   #g   #H   #R    fidelity      alpha      |grad|    ');
        report_optim(param,'--------------------------------------------------------------');
    elseif report_struct(2)==1
        report_optim(param,'========================================================================================');
        report_optim(param,'Iter  #f   #g   #H   #R    fidelity      penalise      total       alpha      |grad|    ');
        report_optim(param,'----------------------------------------------------------------------------------------');
    elseif report_struct(1)==1
        report_optim(param,'==============================================================================');
        report_optim(param,'Iter  #f   #g   #H   #R    fidelity      check           alpha      |grad|    ');
        report_optim(param,'------------------------------------------------------------------------------');
    end
else
    if sum(report_struct)==3
        report_optim(param,'================================================================================================================');
        report_optim(param,'Iter  #f   #g   #H   #R   #T   #O    fidelity      check         penalise      total       alpha      |grad|    ');
        report_optim(param,'----------------------------------------------------------------------------------------------------------------');
    elseif sum(report_struct)==1
        report_optim(param,'========================================================================');
        report_optim(param,'Iter  #f   #g   #H   #R   #T   #O    fidelity      alpha      |grad|    ');
        report_optim(param,'------------------------------------------------------------------------');
    elseif report_struct(2)==1
        report_optim(param,'==================================================================================================');
        report_optim(param,'Iter  #f   #g   #H   #R   #T   #O    fidelity      penalise      total       alpha      |grad|    ');
        report_optim(param,'--------------------------------------------------------------------------------------------------');
    elseif report_struct(1)==1
        report_optim(param,'========================================================================================');
        report_optim(param,'Iter  #f   #g   #H   #R   #T   #O    fidelity      check           alpha      |grad|    ');
        report_optim(param,'----------------------------------------------------------------------------------------');
    end
end

end

% Footer printing function
function footer(param,exitflag,data,report_struct)
if report_struct(3)==0
    if sum(report_struct)==2,       report_optim(param,'------------------------------------------------------------------------------------------------------');
    elseif sum(report_struct)==0,   report_optim(param,'--------------------------------------------------------------');
    elseif report_struct(2)==1,     report_optim(param,'----------------------------------------------------------------------------------------');
    elseif report_struct(1)==1,     report_optim(param,'------------------------------------------------------------------------------');
    end
else
    if sum(report_struct)==3,       report_optim(param,'----------------------------------------------------------------------------------------------------------------');
    elseif sum(report_struct)==1,   report_optim(param,'------------------------------------------------------------------------');
    elseif report_struct(2)==1,     report_optim(param,'--------------------------------------------------------------------------------------------------');
    elseif report_struct(1)==1,     report_optim(param,'----------------------------------------------------------------------------------------');
    end
end

switch(param.method)
    case 'lbfgs',  data.algorithm='LBFGS quasi-Newton method';
    case 'newton-raphson', data.algorithm='Newton-Raphson method';
    case 'gauss-newton', data.algorithm='Gauss-Newton method';
end
switch exitflag
    case  1, message='norm(gradient,2) < tol_gfx';
    case  2, message='norm(step,1) < tol_x';
    case  3, message='fx > tol_f';
    case  0, message='number of iterations exceeded';
    case -2, message='line search found no minimum';
end
report_optim(param,['    Algorithm Used     : ' data.algorithm]);
report_optim(param,['    Exit message       : ' message]);
report_optim(param,['    Iterations         : ' int2str(data.count.iter)]);
report_optim(param,['    Function Count     : ' int2str(data.count.fx)]);
report_optim(param,['    Gradient Count     : ' int2str(data.count.gfx)]);
report_optim(param,['    Hessian Count      : ' int2str(data.count.hfx)]);

if report_struct(3)==0
    if sum(report_struct)==2,       report_optim(param,'======================================================================================================');
    elseif sum(report_struct)==0,   report_optim(param,'==============================================================');
    elseif report_struct(2)==1,     report_optim(param,'========================================================================================');
    elseif report_struct(1)==1,     report_optim(param,'==============================================================================');
    end
else
    if sum(report_struct)==3,       report_optim(param,'================================================================================================================');
    elseif sum(report_struct)==1,   report_optim(param,'========================================================================');
    elseif report_struct(2)==1,     report_optim(param,'==================================================================================================');
    elseif report_struct(1)==1,     report_optim(param,'========================================================================================');
    end
end
end

% Iteration report function
function itrep(param,fx,g,alpha,data,report_struct)

if report_struct(3)==0
    if sum(report_struct)==2
        
        % Performance figures
        fid=data.fx_sep_pen(1);
        pens=sum(-data.fx_sep_pen(2:end));
        chk=data.fx_chk(1);
        
        % Print iteration data
        report_optim(param,[pad(num2str(data.count.iter,'%4.0f'),6),...
            pad(num2str(data.count.fx,'%4.0f'),5),...
            pad(num2str(data.count.gfx,'%4.0f'),5),...
            pad(num2str(data.count.hfx,'%4.0f'),5),...
            pad(num2str(data.count.rfo,'%4.0f'),5),...
            pad(num2str(fid,'%+9.8f'),11),'  '...
            pad(num2str(chk,'(%+11.8f)'),11),'  '...
            pad(num2str(pens,'%+9.6f'),11),'   '...
            pad(num2str(fx,'%+9.6f'),11),'  '...
            pad(num2str(alpha,'%7.2e'),9),'  '...
            pad(num2str(norm(g,2),'%10.4e'),10)],func2str(param.optimcon_fun));
    elseif sum(report_struct)==0
        
        % Performance figures
        fid=data.fx_sep_pen(1);
        
        % Print iteration data
        report_optim(param,[pad(num2str(data.count.iter,'%4.0f'),6),...
            pad(num2str(data.count.fx,'%4.0f'),5),...
            pad(num2str(data.count.gfx,'%4.0f'),5),...
            pad(num2str(data.count.hfx,'%4.0f'),5),...
            pad(num2str(data.count.rfo,'%4.0f'),5),...
            pad(num2str(fid,'%+9.8f'),11),'    '...
            pad(num2str(alpha,'%7.2e'),9),'  '...
            pad(num2str(norm(g,2),'%10.4e'),10)],func2str(param.optimcon_fun));
        
    elseif report_struct(2)==1
        
        % Performance figures
        fid=data.fx_sep_pen(1);
        pens=sum(-data.fx_sep_pen(2:end));
        
        % Print iteration data
        report_optim(param,[pad(num2str(data.count.iter,'%4.0f'),6),...
            pad(num2str(data.count.fx,'%4.0f'),5),...
            pad(num2str(data.count.gfx,'%4.0f'),5),...
            pad(num2str(data.count.hfx,'%4.0f'),5),...
            pad(num2str(data.count.rfo,'%4.0f'),5),...
            pad(num2str(fid,'%+9.8f'),11),'   '...
            pad(num2str(pens,'%+9.6f'),11),'   '...
            pad(num2str(fx,'%+9.6f'),11),'  '...
            pad(num2str(alpha,'%7.2e'),9),'  '...
            pad(num2str(norm(g,2),'%10.4e'),10)],func2str(param.optimcon_fun));
        
    elseif report_struct(1)==1
        
        % Performance figures
        fid=data.fx_sep_pen(1);
        chk=data.fx_chk(1);
        
        % Print iteration data
        report_optim(param,[pad(num2str(data.count.iter,'%4.0f'),6),...
            pad(num2str(data.count.fx,'%4.0f'),5),...
            pad(num2str(data.count.gfx,'%4.0f'),5),...
            pad(num2str(data.count.hfx,'%4.0f'),5),...
            pad(num2str(data.count.rfo,'%4.0f'),5),...
            pad(num2str(fid,'%+9.8f'),11),'   '...
            pad(num2str(chk,'(%+11.8f)'),11),'    '...
            pad(num2str(alpha,'%7.2e'),9),'  '...
            pad(num2str(norm(g,2),'%10.4e'),10)],func2str(param.optimcon_fun));
        
    end
else
    if sum(report_struct)==3
        
        % Performance figures
        fid=data.fx_sep_pen(1);
        pens=sum(-data.fx_sep_pen(2:end));
        chk=data.fx_chk(1);
        
        % Print iteration data
        % HACK this doesn't consider the whole ensemble for trotter number
        report_optim(param,[pad(num2str(data.count.iter,'%4.0f'),6),...
            pad(num2str(data.count.fx,'%4.0f'),5),...
            pad(num2str(data.count.gfx,'%4.0f'),5),...
            pad(num2str(data.count.hfx,'%4.0f'),5),...
            pad(num2str(data.count.rfo,'%4.0f'),5),...
            pad(num2str(param.drift_sys{1}.trotter_number,'%4.0f'),5),...
            pad(num2str(param.drift_sys{1}.split_order,'%4.0f'),5),...
            pad(num2str(fid,'%+9.8f'),11),'  '...
            pad(num2str(chk,'(%+11.8f)'),11),'  '...
            pad(num2str(pens,'%+9.6f'),11),'   '...
            pad(num2str(fx,'%+9.6f'),11),'  '...
            pad(num2str(alpha,'%7.2e'),9),'  '...
            pad(num2str(norm(g,2),'%10.4e'),10)],func2str(param.optimcon_fun));
    elseif sum(report_struct)==1
        
        % Performance figures
        fid=data.fx_sep_pen(1);
        
        % Print iteration data
        % HACK this doesn't consider the whole ensemble for trotter number
        report_optim(param,[pad(num2str(data.count.iter,'%4.0f'),6),...
            pad(num2str(data.count.fx,'%4.0f'),5),...
            pad(num2str(data.count.gfx,'%4.0f'),5),...
            pad(num2str(data.count.hfx,'%4.0f'),5),...
            pad(num2str(data.count.rfo,'%4.0f'),5),...
            pad(num2str(param.drift_sys{1}.trotter_number,'%4.0f'),5),...
            pad(num2str(param.drift_sys{1}.split_order,'%4.0f'),5),...
            pad(num2str(fid,'%+9.8f'),11),'    '...
            pad(num2str(alpha,'%7.2e'),9),'  '...
            pad(num2str(norm(g,2),'%10.4e'),10)],func2str(param.optimcon_fun));
        
    elseif report_struct(2)==1
        
        % Performance figures
        fid=data.fx_sep_pen(1);
        pens=sum(-data.fx_sep_pen(2:end));
        
        % Print iteration data
        % HACK this doesn't consider the whole ensemble for trotter number
        report_optim(param,[pad(num2str(data.count.iter,'%4.0f'),6),...
            pad(num2str(data.count.fx,'%4.0f'),5),...
            pad(num2str(data.count.gfx,'%4.0f'),5),...
            pad(num2str(data.count.hfx,'%4.0f'),5),...
            pad(num2str(data.count.rfo,'%4.0f'),5),...
            pad(num2str(param.drift_sys{1}.trotter_number,'%4.0f'),5),...
            pad(num2str(param.drift_sys{1}.split_order,'%4.0f'),5),...
            pad(num2str(fid,'%+9.8f'),11),'   '...
            pad(num2str(pens,'%+9.6f'),11),'   '...
            pad(num2str(fx,'%+9.6f'),11),'  '...
            pad(num2str(alpha,'%7.2e'),9),'  '...
            pad(num2str(norm(g,2),'%10.4e'),10)],func2str(param.optimcon_fun));
        
    elseif report_struct(1)==1
        
        % Performance figures
        fid=data.fx_sep_pen(1);
        chk=data.fx_chk(1);
        
        % Print iteration data
        % HACK this doesn't consider the whole ensemble for trotter number
        report_optim(param,[pad(num2str(data.count.iter,'%4.0f'),6),...
            pad(num2str(data.count.fx,'%4.0f'),5),...
            pad(num2str(data.count.gfx,'%4.0f'),5),...
            pad(num2str(data.count.hfx,'%4.0f'),5),...
            pad(num2str(data.count.rfo,'%4.0f'),5),...
            pad(num2str(param.drift_sys{1}.trotter_number,'%4.0f'),5),...
            pad(num2str(param.drift_sys{1}.split_order,'%4.0f'),5),...
            pad(num2str(fid,'%+9.8f'),11),'   '...
            pad(num2str(chk,'(%+11.8f)'),11),'    '...
            pad(num2str(alpha,'%7.2e'),9),'  '...
            pad(num2str(norm(g,2),'%10.4e'),10)],func2str(param.optimcon_fun));
        
    end
end
end

function report_optim(param,report_string,optimcon_fun)

% List call stack exceptions
not_useful={'make_general_channel/channel_general','parallel_function','header','footer','itrep','report_optim'};

% Compose the prefix
call_stack=dbstack;
for n=1:numel(call_stack)
    if ismember(call_stack(n).name,not_useful), call_stack(n).name='';
    elseif ismember(call_stack(n).name,{'remoteParallelFunction','remoteBlockExecution'}), call_stack(n).name='worker_node > ';
    else, call_stack(n).name=[call_stack(n).name ' > '];
    end
end

prefix_string=[call_stack(end:-1:2).name]; prefix_string=prefix_string(1:end-3);
if nargin==3, prefix_string=[prefix_string '@' optimcon_fun]; end

% Fix empty prefixes
if isempty(prefix_string), prefix_string=' '; end

% Roll the prefix
if numel(prefix_string)<50, prefix_string=pad(prefix_string,50);
else, prefix_string=['...' prefix_string((end-46):end)];
end

% Send the report string to the output, ignoring impossible writes
try fprintf(param.output,'%s\n',['[' prefix_string ' ]  ' report_string]); end %#ok<TRYNC>

end

function [H,data]=hessianreg(param,H,g,data)

% Set shorthands
alpha=param.reg_alpha;

% Check if RFO is needed
[~,p]=chol(H); pos_def=~logical(p);
if pos_def&&(cond(H,2)<param.reg_max_cond), return; end

% Start RFO iteration loop
for ind=1:param.reg_max_iter
    
    % Make auxiliary Hessian
    H=[(alpha^2)*H   alpha*g;
        alpha*g'     0];
    
    % Calculate eigenvalue shift
    sigma=min([0 min(eig(H,'vector'))]);
   
    % Apply eigenvalue shift
    H=H-sigma*speye(size(H));
    
    % Truncate and scale back
    H=H(1:(end-1),(1:end-1)); H=H/(alpha^2);
    
    % Update alpha
    alpha=alpha*param.reg_phi;
    
    % Increment the counter
    data.count.rfo=data.count.rfo+1;
    
    % Break if condition number reached
    if cond(H,2)<param.reg_max_cond, break; end
    
end

% Clean up the result
H=real(H+H')/2;

end

function direction=search_direction(method,dx_hist,dg_hist,g)

if strcmp(method,'lbfgs')

    % Initialize variables
    N=size(dx_hist,2);
    alpha=zeros(1,N);
    p=zeros(1,N);

    % Loop over history
    for n=1:N
        p(n)=1/(dg_hist(:,n)'*dx_hist(:,n));
        alpha(n)=p(n)*dx_hist(:,n)'*g;
        g=g-alpha(n)*dg_hist(:,n);
    end

    % Scaling of initial Hessian (identity matrix)
    p_k=dg_hist(:,1)'*dx_hist(:,1)/sum(dg_hist(:,1).^2);

    % Make r = - Hessian * gradient
    direction=p_k*g;
    for n=N:-1:1
        b=p(n)*dg_hist(:,n)'*direction;
        direction=direction+dx_hist(:,n)*(alpha(n)-b);
    end
    % code above for minimisation
    direction=-direction;

elseif strcmp(method,'steep')
    direction=-g;
end

end

% Consistency enforcement
function grumble(optim_param,cost_function,guess)
if ~isa(cost_function,'function_handle'), error('cost_function must be a function handle.'); end
if all(guess==0), warning('Bad things happen when the initial guess is zeros.'); end
end
