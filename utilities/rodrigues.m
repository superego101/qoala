% Function to calculate propagator elements using the Euler-Rodrigues
% formula, instead of an expensive matrix exponential. The function also
% calculates derivatives and second derivatives, with respect to control 
% operatos, instead of using the auxiliary matrix formalism.
%
% Author:       David L. Goodwin 
% Dev. period:  01/2020--04/2022
% Updated:      13/04/22
%
% Contact:
%   david.goodwin@chem.ox.ac.uk
%   david.goodwin@partner.kit.edu

function [P,dP,d2P]=rodrigues(state_space,basis,rot_ops,ops_axis,rot_axis,dt,method,spin_control)

% check number of outputs required
n_outputs=nargout;

% shorthand for Cartesian controls
x=rot_axis(:,1); y=rot_axis(:,2); z=rot_axis(:,3);

% find the amplitude
r = sqrt(x.^2+y.^2+z.^2);

% shorthand for rotation angle
delta=dt.*r;

% components of the quaternion from the xyz rotation axis and angle
rot_axis=[x, y, z];

% Normalize the axis vector
rot_axis=rot_axis./r;

% remove NaN from zeros
rot_axis(isnan(rot_axis))=0;

% shorthand complex numbers
sind=sin(delta/2);
alpha=cos(delta/2) - 1i*rot_axis(:,3).*sind;
beta= (rot_axis(:,2) - 1i*rot_axis(:,1)).*sind;

% construct the propagators using Rodrigues' formula
if strcmp(state_space,'liouville') && strcmp(basis,'zeeman')
        
        % shortcuts
        alpha2=alpha.^2;
        beta2=beta.^2;
        alphabeta=(alpha.*beta);
        alphaconjbeta=(conj(alpha).*beta);
        alphabetaconj=(alpha.*conj(beta));
        alphaconjalpha=(conj(alpha).*alpha);
        betaconjbeta=(conj(beta).*beta);

        P=[+alphaconjalpha, +alphaconjbeta,  +alphabetaconj,   +betaconjbeta,...
            -conj(alphabeta),  +conj(alpha2),     -conj(beta2),       +conj(alphabeta),...
            -alphabeta,        -beta2,            +alpha2,            +alphabeta,...
            +betaconjbeta,   -alphaconjbeta,  -alphabetaconj,   +alphaconjalpha];
        
elseif strcmp(state_space,'liouville') && strcmp(basis,'sphten')
        
        % shortcuts
        alpha2=alpha.^2;
        beta2=beta.^2;
        alphabeta=sqrt(2)*(alpha.*beta);
        alphaconjbeta=sqrt(2)*(conj(alpha).*beta);
        
        % (Wigner matrix)
        P=[alpha2,                   alphabeta,                  beta2,...
            -conj(alphaconjbeta),  abs(alpha2)-abs(beta2),  alphaconjbeta,...
            conj(beta2),              -conj(alphabeta),           conj(alpha2)];
        
elseif strcmp(state_space,'hilbert') && strcmp(basis,'zeeman')
        
        P=[+alpha,      +beta,...
            -conj(beta), +conj(alpha)];
        
end

if n_outputs>1
    
    % cleanup zeros
    r(r==0)=eps;
    
    % precalculate trig functions
    cos_norm=cos(delta);
    sin_norm=sin(delta);
    
    % zero vector
    o=zeros(length(z),1);
    
    % shorthand squares
    xy=x.*y; xz=x.*z; yz=y.*z;
    x2=x.^2; y2=y.^2; z2=z.^2; r2=r.^2; %dtr2=dt.*r2;    
    
    T = [    dt,      o,      o,...
              o,     dt,      o];%;...
             %o;      o;     dt];
    
    S = [     o,     +z,     -y,...
             -z,      o,     +x];%;...
            %+y;     -x;      o];
    
    S2= [-y2-z2,    +xy,    +xz,...
            +xy, -x2-z2,    +yz];%;...
           %+xz;    +yz; -x2-y2];
    
%     T=[]; S=[]; S2=[]; % initialise
%     % should be unique(ops_axis), otherwise there is repeated calculation
%     for k=1:sum(spin_control())
%     for k=1:numel(ops_axis)
%         if strcmp(ops_axis{k},'x')
%             T =[T,      dt,      o,      o];
%             S =[S,       o,     +z,     -y];
%             S2=[S2, -y2-z2,    +xy,    +xz];
%         elseif strcmp(ops_axis{k},'y')
%             T =[T,       o,     dt,      o];
%             S =[S,      -z,      o,     +x];
%             S2=[S2,    +xy, -x2-z2,    +yz];
%         elseif strcmp(ops_axis{k},'z')
%             T =[T,       o,      o,     dt];
%             S =[S,      +y,     -x,      o];
%             S2=[S2,    +xz,    +yz, -x2-y2];
%         end
%     end
    
    c1 = (cos_norm-1)./r2;
    c2 = (dt-sin_norm./r)./r2;
    
    D = (T +c1.*S + c2.*S2);
    
    % find the number of spins, controls (x/y), and time steps
    nspins=size(rot_ops,1);
    kctrls=2*size(rot_ops,2)/3;
    nsteps=length(dt)/nspins;
    
    % initialise the control derivative cell
    dP=cell(1,kctrls);
    
    % initialise the spin counter
    spin_no=1:nspins;
    
    for k=1:kctrls/2
        
         % decide if control affects more than one spin operator
         if ~isempty(spin_control) && sum(spin_control(:,2*k-1))>1
            
            if strcmp(method,'full')
                % initialisations
                sigma=[]; D_stack=[]; ind=[]; counter=0;
                
                % loop over all spins
                for s=spin_no
                    
                    % project out single spin parts of control operator
                    if spin_control(s,2*k-1)
                        
                        % increment the counters
                        counter=counter+1;
                        
                        % pull the (single-spin) Cartesian control operators
                        Lx=rot_ops{s,3*k-2};
                        Ly=rot_ops{s,3*k-1};
                        Lz=rot_ops{s,3*k};
                        
                        % stack sigma matrix, derivative matrix, and indices
                        sigma   = [sigma; -1i*[Lx(:), Ly(:), Lz(:)].'];       %#ok<AGROW>
                        D_stack = [D_stack, D((s-1)*nsteps +1:s*nsteps,:)]; %#ok<AGROW>
                        ind     = [ind 1+(counter-1)*6:3+(counter-1)*6];    %#ok<AGROW>
                        
                        % remove the spin from the list
                        spin_no(spin_no==s)=[];
                        
                    end
                    
                end
                
                % calculate the propagator derivatives
                dP{2*k-1} = sparse(D_stack(:,ind))*sigma;
                dP{2*k}   = sparse(D_stack(:,ind+3))*sigma;
                %dP{3} = sigma*sparse(D_stack(ind+6,:));
                
%             elseif strcmp(method,'sparse')
%                 
%                 % initialisations
%                 sigma=[]; D_stack=[]; ind=[]; counter=0;
%                 
%                 % loop over all spins
%                 for s=spin_no
%                     
%                     % project out single spin parts of control operator
%                     if spin_control(s,2*k-1)
%                         
%                         % increment the counters
%                         counter=counter+1;
%                         
%                         % pull the (single-spin) Cartesian control operators
%                         Lx=rot_ops{s,3*k-2};
%                         Ly=rot_ops{s,3*k-1};
%                         Lz=rot_ops{s,3*k};
%                         
%                         % stack sigma matrix, derivative matrix, and indices
%                         sigma   = -1i*[Lx(:), Ly(:), Lz(:)];       %#ok<AGROW>
%                         D_stack = D(:,(s-1)*nsteps +1:s*nsteps); %#ok<AGROW>
%                         
%                         % calculate the propagator derivatives
%                         dP{2*k-1}(:,(counter-1)*nsteps +1:counter*nsteps) = sigma*sparse(D_stack([1,2,3],:));
%                         dP{2*k}(:,(counter-1)*nsteps +1:counter*nsteps)   = sigma*sparse(D_stack([1,2,3]+3,:));
%                         %dP{3} = sigma*sparse(D_stack(ind+6,:));
%                         
%                         % remove the spin from the list
%                         spin_no(spin_no==s)=[];
%                         
%                     end
%                     
%                 end
%                 
             end
             
         else
            
            % pull the (single-spin) Cartesian control operators
            Lx=rot_ops{spin_no(1),3*k-2};
            Ly=rot_ops{spin_no(1),3*k-1};
            Lz=rot_ops{spin_no(1),3*k};
            
            % construct 'sigma' transform matrix
            sigma = sparse(-1i*[Lx(:) Ly(:) Lz(:)].');
            
            % calculate the propagator derivatives
            dP{2*k-1} = sparse(D((k-1)*nsteps +1:k*nsteps,[1,2,3]))*sigma; % SPARSE HACK!!
            dP{2*k}   = sparse(D((k-1)*nsteps +1:k*nsteps,[4,5,6]))*sigma;
            %dP{3} = sigma*D([7,8,9],(k-1)*nsteps +1:k*nsteps);
            
            % remove the spin from the list
            spin_no(1)=[];
         end
         
    end
    
    if n_outputs>2
        
        % HESSIAN CODE NOT YET GENERAL FOR MULTI-SPIN
        dc1 = ((2-2*cos_norm)./(delta) - sin_norm)./(r2.*r);
        dc2 = (-2-cos_norm + (3*sin_norm)./(delta))./(r2.*r2);
        
        %SaI = [     o;      o;    +c1;...
        %            o;      o;    -c1];%;...
        %          % ?;      ?;      ?];
        
        SaSx= [     o;     +y;     +z;...
                   +y;   -2*x;      o];%;...
                  % ?;      ?;      ?];
        
        SaSy= [  -2*y;     +x;      o;...
                   +x;      o;     +z];%;...
                  % ?;      ?;      ?];
        
        %SaSz= [     ?;      ?;      ?;...
        %            ?;      ?;      ?;...
        %            ?;      ?;      ?];
        
        dDx = dt.*(dc1.*x.*S1 + dc2.*x.*S2) + (c2.*SaSx + c1.*[0;0;0;0;0;1])./r2;
        dDy = dt.*(dc1.*y.*S1 + dc2.*y.*S2) + (c2.*SaSy - c1.*[0;0;1;0;0;0])./r2;
        %dDz = dt.*(dc1.*z.*S1 + dc2.*z.*S2) + (c2.*SaSz - c1.*[?;?;?;?;?;?])./r2;
        
        %dDx(6,:) = dDx(6,:) + c1./r2;
        %dDy(3,:) = dDy(3,:) - c1./r2;
        
        d2P{1,1} = sigma*(dDx([1,2,3],:));
        d2P{1,2} = sigma*(dDy([1,2,3],:));
        %d2P{1,3} = sigma*(dDx([1,2,3],:));
        d2P{2,1} = sigma*(dDx([4,5,6],:));
        d2P{2,2} = sigma*(dDy([4,5,6],:));
        %d2P{2,3} = sigma*(dDz([4,5,6],:));
        %d2P{3,1} = sigma*(dDx([7,8,9],:));
        %d2P{3,2} = sigma*(dDy([7,8,9],:));
        %d2P{3,3} = sigma*(dDx([7,8,9],:));
        
        if n_outputs>3
            
            error('incorrect number of outputs')
        end
    end
end


end
