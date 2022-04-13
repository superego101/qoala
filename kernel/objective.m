% Calls and collect the correct amount of outputs from an objective 
% function - used by optimisation routines. Syntax:
%
%  [data,fx,grad,hess]=objective(x,objfun_handle,data,ctrl_param)
%
% Parameters:
%
%    x                  - objective function argument 
%    
%    objfun_handle      - handle to the objective function
%    
%    data               - data structure inherited from 
%                         fminnewton.m
%
% Outputs:
%
%    data               - modified data structure with 
%                         diagnostics from the objective
%    
%    fx                 - objective function value at x
%    
%    grad               - gradient of the objective function 
%                         at x
%    
%    hess               - Hessian of the objective function 
%                         at x
%
% Author:       David L. Goodwin 
% Dev. period:  01/2020--04/2022
% Updated:      13/04/22
%
% Contact:
%   david.goodwin@chem.ox.ac.uk
%   david.goodwin@partner.kit.edu
%
% NOTE: This function is based on the Spinach (Goodwin, Kuprov) function 
%       objeval(), and should be referenced accordingly.

function [data,ctrl_param,fx,grad,hess]=objective(x,objfun_handle,data,ctrl_param)

% Check consistency
grumble(objfun_handle,x)

% start the objective function timer
obj_time=tic;

% Switch between function/gradient/Hessian calls
switch nargout
    
    case 3
        
        % Run the objective function for the fidelity
        [traj_data,fidelity]=feval(objfun_handle,reshape(x,data.x_shape),ctrl_param);
        
        % Increment counters
        data.count.fx=data.count.fx+1;
        
        % Assign the data structure
        data.fx_sep_pen=fidelity;
        data.traj_data=traj_data;
        
        % Compute total error functional
        fx=fidelity(1)-sum(fidelity(2:end));
    
    case 4
        
        
        % Run the objective function for the fidelity and gradient
        [traj_data,fidelity,grad]=feval(objfun_handle,reshape(x,data.x_shape),ctrl_param);
        
        % Increment counters
        data.count.fx=data.count.fx+1;
        data.count.gfx=data.count.gfx+1;
        
        % Assign the data structure
        data.fx_sep_pen=fidelity;
        data.traj_data=traj_data;
        
        % Compute total error functional
        fx=fidelity(1)-sum(fidelity(2:end));
        
        % Compute the total gradient
        grad=grad(:,:,1)-sum(grad(:,:,2:end),3); grad=grad(:);
    
    case 5
        
        % Run the objective function for the fidelity, gradient, and Hessian
        [traj_data,fidelity,grad,hess]=feval(objfun_handle,reshape(x,data.x_shape),ctrl_param);
        
        % Increment counters
        data.count.fx=data.count.fx+1;
        data.count.gfx=data.count.gfx+1;
        data.count.hfx=data.count.hfx+1;
        
        % Assign the data structure
        data.fx_sep_pen=fidelity;
        data.traj_data=traj_data;
        
        % Compute total error functional
        fx=fidelity(1)-sum(fidelity(2:end));
        
        % Compute the total gradient
        grad=grad(:,:,1)-sum(grad(:,:,2:end),3); grad=grad(:);
        
        % Compute the total Hessian
        hess=hess(:,:,1)-sum(hess(:,:,2:end),3);
        
    otherwise
        
        % Complain and bomb out
        error('must have at least two output arguments.');
        
end

% find the current iteration
iter=find(~isnan(data.timer(:,1)),1,'last');
if isnan(data.timer(iter,3))
    
    % set as the timer
    if iter==1
        data.timer(iter,3)=toc(obj_time);
    else
        data.timer(iter,3)=data.timer(iter-1,3)+toc(obj_time);
    end
else
    % add to the timer
    data.timer(iter,3)=data.timer(iter,3)+toc(obj_time);
end
if isfield(data.traj_data,'timer')
    if isnan(data.timer(iter,4))

        % set as the timer
        if iter==1
            data.timer(iter,4)=data.traj_data.timer;
        else
            data.timer(iter,4)=data.timer(iter-1,4)+data.traj_data.timer;
        end
    else
        % add to the timer
        data.timer(iter,4)=data.timer(iter,4)+data.traj_data.timer;
    end
end

% update the spin_system structure
if isfield(traj_data,'drift_sys')
    for n=1:numel(traj_data.drift_sys)
        ctrl_param.drift_sys{n}=traj_data.drift_sys{n};
    end
end

if isfield(ctrl_param,'fidelity_chk')
    ctrl_param_chk=ctrl_param;
    ctrl_param_chk.optimcon_fun=ctrl_param.fidelity_chk;
    [~,data.fx_chk]=feval(objfun_handle,reshape(x,data.x_shape),ctrl_param_chk);
else
    data.fx_chk=fidelity(1);
end

if isfield(traj_data,'f_adapt_count')
    data.count.fx=data.count.fx+traj_data.f_adapt_count;
    traj_data.f_adapt_count=0;
end

if isfield(traj_data,'grad_norm')
    ctrl_param.grad_norm=traj_data.grad_norm;
end

end

% Consistency enforcement
function grumble(objfun_handle,x)
if ~isa(objfun_handle,'function_handle')
    error('objfun_handle must be a function handle.')
end
if isempty(x)||(~isnumeric(x))||(~isreal(x))
    error('x must be a vector of real numbers')
end
end
