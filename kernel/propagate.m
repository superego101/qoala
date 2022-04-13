% Calculates exponential propagator expm(-1i*L*timestep). Syntax:
%
%           P=propagator(state_space,L,timestep,method,nonzeros_tol,rho)
%
% Inputs:
%
%    state_space-  either 'liouville' or 'hilbert'
%
%    L          -  Hamiltonian or Liouvillian matrix
%
%    timestep   -  propagation time step
%
%    method     -  expm method, options: 'pade' for Matlabs Pade
%                  approximation, 'krylov' for krylov propagation (requires
%                  rho as input), 'taylor' for a truncated Taylor series
%                  with scaling and squaring.
%
%    nonzero_tol-  remove elements below this value
%
%    rho        -  state to propagate from (optional)
%
% Output:
%
%    P          -  propagator matrix or propagated state, if rho is
%                  provided.
%
% NOTE: 'krylov' uses a form of the spinach function step.m
% NOTE: 'taylor' uses a form of the spinach function propagator.m
%
% Author:       David L. Goodwin 
% Dev. period:  01/2020--04/2022
% Updated:      13/04/22
%
% Contact:
%   david.goodwin@chem.ox.ac.uk
%   david.goodwin@partner.kit.edu
%
% NOTE: This function is based on the Spinach function propagator(), and
%       should be referenced accordingly.

function P=propagate(state_space,L,timestep,method,nonzero_tol,rho)

% Set a shorthand for -i*L*dt
A=-1i*L*timestep;

% return state or propagator
if nargin>5, state_prop=true; else, state_prop=false; end

if strcmp(method,'pade')
    
    % matlab's Pade approximant
    P=expm(full(A)); % ensure full matrix

    if state_prop && strcmp(state_space,'liouville') 
        P=P*rho; 
    elseif state_prop && strcmp(state_space,'hilbert') 
        P=P*rho*P';
    end

    % re-assert sparsity
    if issparse(L), P=sparse(P); end

elseif strcmp(method,'taylor')

    % Estimate the norm
    mat_norm=norm(full(A),2); % Spinach' cheap_norm() better
    
    % Determine scaling and squaring parameters
    n_squarings=max([0 ceil(log2(mat_norm))]); scaling_factor=2^n_squarings;

    % Scale and clean up the matrix
    if scaling_factor>1, A=(1/scaling_factor)*A; end
    A=nonzero_tol*round((1/nonzero_tol)*(A));

    % Run Taylor series procedure on the CPU
    P=speye(size(A)); next_term=P; n=1;
    
    while nnz(next_term)>0

        % Compute the next term
        if issparse(A)
            next_term=(1/n)*A*next_term;
        else
            next_term=(1/n)*next_term*A;
        end

        % Eliminate small elements
        next_term=nonzero_tol*round((1/nonzero_tol)*(next_term));

        % Add to the total and increment the counter
        P=P+next_term; n=n+1;

    end

    % Reclaim memory
    clear('A','next_term');

    % Clean up the result
    P=nonzero_tol*round((1/nonzero_tol)*(P));

    % Run the squaring stage
    if n_squarings>0

        % Run serial CPU squaring
        for n=1:n_squarings

            % Square the propagator and clean
            P=nonzero_tol*round((1/nonzero_tol)*(P*P));

        end

    end

    if state_prop && strcmp(state_space,'liouville') 
        P=P*rho; 
    elseif state_prop && strcmp(state_space,'hilbert') 
        P=P*rho*P'; 
    end

elseif strcmp(method,'krylov')

    % Make sure rho is full
    if issparse(rho), rho=full(rho); end

    % Determine the number of time steps
    norm_mat=norm(full(L),2)*abs(timestep);  % Spinach' cheap_norm() better
    nsteps=ceil(norm_mat/2);

    % Estimate the scaling coefficient
    scaling=max([1 norm(rho,2)]);  % Spinach' cheap_norm() better

    % Scale the vector
    rho=rho/scaling;

    % Run the Krylov procedure
    for n=1:nsteps
        next_term=rho; k=1;
        while nnz(abs(next_term)>eps)>0
            next_term=(-1i*timestep/(k*nsteps))*(L*next_term);
            rho=rho+next_term; k=k+1;
        end
    end

    % Scale the vector back
    P=scaling*rho;

else
    error(['expm method ''' method ''' not coded'])
end

end
