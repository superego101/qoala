% Penalty terms for the Optimal Control module. Returns the penalty
% function and its gradient for the waveform.
%
% Author:       David L. Goodwin 
% Dev. period:  01/2020--04/2022
% Updated:      13/04/22
%
% Contact:
%   david.goodwin@chem.ox.ac.uk
%   david.goodwin@partner.kit.edu
%
% NOTE: This function is rigidly based on the Spinach (Goodwin, Kuprov) 
%       function objeval(), and should be referenced accordingly.

function [pen_term,pen_grad,pen_hess]=penalty_fun(param,wf,type,fb,cb,weight)

% Preallocate the results
if nargout>0, pen_term=0; end
if nargout>1, pen_grad=zeros(size(wf)); end
if nargout>2, pen_hess=zeros(numel(wf),numel(wf)); end

% Decide the penalty type
switch type
    
    case 'none'
        
        % Nothing to do
        
    case 'NS'
        
        % Compute the penalty
        if nargout>0
            pen_term=sum(sum(wf.^2));
            pen_term=pen_term/size(wf,2);
        end
        
        % Compute the gradient
        if nargout>1
            pen_grad=2*wf;
            pen_grad=pen_grad/size(wf,2);
        end
        
        % Compute the Hessian
        if nargout>2
            pen_hess=2*speye(numel(wf));
            pen_hess=pen_hess/size(wf,2);
        end
        
    case 'DNS'
        
        % Five-point second derivative matrix
        D=fdmat(size(wf,2),5,2,'wall');
        
        % Compute the penalty
        if nargout>0
            dwf=wf*D';
            pen_term=norm(dwf,2)^2;
            pen_term=pen_term/size(dwf,2);
        end
        
        % Compute the gradient
        if nargout>1
            pen_grad=2*dwf*D;
            pen_grad=pen_grad/size(dwf,2);
        end
        
        % Compute the Hessian
        if nargout>2
            pen_hess=2*kron(D'*D,speye(size(wf,1)));
            pen_hess=pen_hess/size(wf,2);
        end
        
    case 'SNS'
        
        % Build ceiling hole inventory
        ch_map=(wf>cb); ch_actual=wf.*ch_map; ch_wanted=cb.*ch_map;
        
        % Build floor hole inventory
        fh_map=(wf<fb); fh_actual=wf.*fh_map; fh_wanted=fb.*fh_map;
        
        % Compute the penalty
        pen_term=pen_term+sum(sum((ch_actual-ch_wanted).^2))+...
            sum(sum((fh_actual-fh_wanted).^2));
        pen_term=pen_term/size(wf,1);
        
        % Compute the gradient
        if nargout>1
            pen_grad=2*(ch_actual-ch_wanted)+...
                2*(fh_actual-fh_wanted);
            pen_grad=pen_grad/size(wf,1);
        end
        
        % Compute the Hessian
        if nargout>2
            pen_hess=2*ch_map+2*fh_map;
            pen_hess=diag(pen_hess(:));
            pen_hess=pen_hess/size(wf,1);
        end
        
    case 'SNSA'
        
        % preallocate amplitude and phase vectors
        amp=zeros(size(wf,1)/2,size(wf,2));
        phi=zeros(size(wf,1)/2,size(wf,2));
        
        % calculate the amplitude and phase of the waveform
        for n=1:size(wf,1)/2
            [amp(n,:),phi(n,:)]=cartesian2polar(wf(2*n-1,:),wf(2*n,:));
        end
        
        % find the amplitude hole inventory
        ch_map=(amp>cb); ch_wanted=cb.*ch_map;
        
        % preallocated new cartesian bounds
        b=zeros(size(wf));
        
        % transform map to cartesian representation, giving new bounds
        for n=1:size(wf,1)/2
            [b(2*n-1,:),b(2*n,:)]=polar2cartesian(ch_wanted(n,:),phi(n,:));
        end
        
        % new cartesian bounds from the amplitude bounds
        cb=b;  cb(cb<=0)=max(max(wf))+1;
        fb=b;  fb(fb>=0)=min(min(wf))-1;
        
        % recursively call the SNS penalty with new bounds
        if nargout==1
            [pen_term]=penalty(wf,'SNS',fb,cb);
        elseif nargout==2
            [pen_term,pen_grad]=penalty(wf,'SNS',fb,cb);
        elseif nargout==3
            [pen_term,pen_grad,pen_hess]=penalty(wf,'SNS',fb,cb);
        end
        
    case 'SNSM'
        
        % transform matrix for gradient
        A=permute(sign(wf),[2,1]);
        A=spdiags(A,-(0:size(wf,2):(size(wf,1)-1)*size(wf,2)),numel(wf),size(wf,2));
        
        % transform the penalty
        wfm=sum(abs(wf),1);
        
        % preallocate amplitude and phase vectors
        % Build ceiling hole inventory
        cb=0.2; fb=0; % hard coded hack
        ch_map=(wfm>cb); ch_actual=wfm.*ch_map; ch_wanted=cb.*ch_map;
        
        % Build floor hole inventory
        fh_map=(wfm<fb); fh_actual=wfm.*fh_map; fh_wanted=fb.*fh_map;
        
        % Compute the penalty
        pen_term=pen_term+sum(sum((ch_actual-ch_wanted).^2))+...
            sum(sum((fh_actual-fh_wanted).^2));
        pen_term=pen_term/size(wfm,2);
        
        % Compute the gradient
        if nargout>1
            pen_grad=2*(ch_actual-ch_wanted)+...
                2*(fh_actual-fh_wanted);
            pen_grad=pen_grad/size(wf,2);
            pen_grad=reshape(pen_grad*A.',[size(wf,2) size(wf,1)]).';
            pen_grad=pen_grad.*wf./sum(wf,1);
        end
        
        % Compute the Hessian
        if nargout>2
            pen_hess=2*ch_map+2*fh_map;
            pen_hess=diag(pen_hess(:));
            pen_hess=pen_hess/size(wf,2);
        end
        
    case 'ADIAB'
        pen_term=0;
        pen_grad=zeros(size(wf));
        pulse_syms=param.pulse_syms;
        
        if (pulse_syms(1)==1 && pulse_syms(2)==1) ||...
                (pulse_syms(1)==-1 && pulse_syms(2)==1) ||...
                (pulse_syms(1)==1 && pulse_syms(2)==-1) ||...
                (pulse_syms(1)==-1 && pulse_syms(2)==-1)
            wf_tmp(1,:)=[wf(1,:) pulse_syms(1)*fliplr(wf(1,:))];
            wf_tmp(2,:)=[wf(2,:) pulse_syms(2)*fliplr(wf(2,:))];
            wf=wf_tmp; clear wf_tmp;
        end
        
        % spread the offsets over the bandwidth with
        offsets=2*pi*param.bandwidth*linspace(-0.5,+0.5,size(wf,2));
        
        for k=1:length(param.pwr_levels)
            
            % pull the Cartesian parts of the effective field
            x=param.pwr_levels(k)*wf(1,:);
            y=param.pwr_levels(k)*wf(2,:);
            z=offsets;
            rot_axis=[x; y; z].';
            
            % Normalize the axis vector
            rot_axis=rot_axis./vecnorm(rot_axis,2,2);
            
            % transform of the effective field to spherical coordinates
            [r_e,theta_e,phi_e] = cartesian2spherical(rot_axis(:,1),rot_axis(:,2),rot_axis(:,3));
            
            % penalise the roughness of the polar angle
            [pen_term_lcl,pen_grad_lcl]=penaltyox(param,(theta_e).','DNS',-1,1);
            pen_term=pen_term+pen_term_lcl;
            
            % transform the gradient back to Cartesian coordinates
            [x,y,z,Dx,Dy,Dz]=spherical2cartesian(r_e.',theta_e.',phi_e.',0*r_e.',pen_grad_lcl,0*phi_e.');
            
            % concatinate
            pen_grad_lcl=[Dx;Dy];
            
            if (pulse_syms(1)==1 && pulse_syms(2)==1) ||...
                    (pulse_syms(1)==-1 && pulse_syms(2)==1) ||...
                    (pulse_syms(1)==1 && pulse_syms(2)==-1) ||...
                    (pulse_syms(1)==-1 && pulse_syms(2)==-1)
                pen_grad_lcl=(pen_grad_lcl(:,1:end/2)+pen_grad_lcl(:,end:-1:end/2+1));
            end
            pen_grad=pen_grad+pen_grad_lcl;
            
        end
        pen_term=pen_term/length(param.pwr_levels);
        pen_grad=pen_grad/length(param.pwr_levels);
        
    otherwise
        
        error('unknown penalty function type.');
        
end

if nargout>0, pen_term=weight*pen_term; end
if nargout>1, pen_grad=weight*pen_grad; end
if nargout>2, pen_hess=weight*pen_hess; end

end
