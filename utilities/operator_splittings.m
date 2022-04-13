% function to output operator splitting propagators and constants for a
% given interaction Hamiltonian, time-step, splitting order, and Trotter
% number.
%
% inputs:
%
%   L1_in   -	interaction Hamiltonian
%   dt      -	time slice
%   odr     -   splitting order
%   trt     -   Trotter number
%
% outputs:
%
%   c_spn   -   splitting constants for the offset Hamiltonians
%   i_spn   -   ordered indices to multiply the offset propagators
%   i_int   -   ordered indices to multiply the interaction propagators
%   P       -   cell array of cell arrays, giving propagators for the
%               interaction Hamiltonians
%
% see:
% Applying splitting methods with complex coefficients to the numerical integration of unitary problems, Blanes 2021
% Splitting methods with complex coefficients for some classes of evolution equations, Blanes
% Splitting and composition methods with embedded error estimators, Blanes 2019
% Splitting methods with complex times for parabolic equations, Castella 2009
%
% Author:       David L. Goodwin 
% Dev. period:  01/2020--04/2022
% Updated:      13/04/22
%
% Contact:
%   david.goodwin@chem.ox.ac.uk
%   david.goodwin@partner.kit.edu

function [c_spn,i_spn,i_int,P]=operator_splittings(param,odr,trt,dt,L1_in)

if odr==0

    i_spn= 1;
    i_int= 1;
    for j=1:trt-1
        i_spn=[i_spn 1]; %#ok<AGROW>
        i_int=[i_int 1]; %#ok<AGROW>
    end
    i_int=[i_int 1];

    c_spn(1)=1;
    c_cpl(1)=0;

elseif odr==1

    i_spn= 1;
    i_int= 1;
    for j=1:trt-1
        i_spn=[i_spn 2]; %#ok<AGROW>
        i_int=[i_int 3]; %#ok<AGROW>
    end
    i_int=[i_int 2];

    c_spn(1)=1;
    c_cpl(1)=0;
    c_cpl(2)=1;

elseif odr==2

    i_spn=1;
    i_int=1;
    for j=1:trt-1
        i_spn=[i_spn 1]; %#ok<AGROW>
        i_int=[i_int 2]; %#ok<AGROW>
    end
    i_int=[i_int 1];

    c_spn(1)=1;
    c_cpl(1)=0.5;

elseif odr==3

    % 1st variant
    i_spn=[1 2 1];
    i_int=[1 2 3];
    for j=1:trt-1
        i_spn=[i_spn 1 2 1]; %#ok<AGROW>
        i_int=[i_int 5 2 3]; %#ok<AGROW>
    end
    i_int=[i_int 4];

    c_spn(1)=+3/10;
    c_spn(2)=+2/5;
    c_cpl(1)=+(13/126) + 1i*((sqrt(59/2))/63);
    c_cpl(2)=+(25/63) - 1i*(5*(sqrt(59/2))/126);
    c_cpl(3)=+(25/63) + 1i*(5*(sqrt(59/2))/126);
    c_cpl(4)=+(13/126) - 1i*((sqrt(59/2))/63);

elseif odr==4

    %variant 1
    i_spn=[1 2 3 3 2 1];
    i_int=[1 2 3 4 3 2];
    for j=1:trt-1
        i_spn=[i_spn 1 2 3 3 2 1]; %#ok<AGROW>
        i_int=[i_int 5 2 3 4 3 2]; %#ok<AGROW>
    end
    i_int=[i_int 1];

    c_spn(1)=+0.209515106613362;
    c_spn(2)=-0.143851773179818;
    c_spn(3)=0.5-(c_spn(1)+c_spn(2));
    c_cpl(1)=+0.0792036964311957;
    c_cpl(2)=+0.353172906049774;
    c_cpl(3)=-0.0420650803577195;
    c_cpl(4)=1-2*(c_cpl(1)+c_cpl(2)+c_cpl(3));

elseif odr==6

    %variant 1
    i_spn=[1 2 3 4 5 5 4 3 2 1];
    i_int=[1 2 3 4 5 6 5 4 3 2];
    for j=1:trt-1
        i_spn=[i_spn 1 2 3 4 5 5 4 3 2 1]; %#ok<AGROW>
        i_int=[i_int 7 2 3 4 5 6 5 4 3 2]; %#ok<AGROW>
    end
    i_int=[i_int 1];

    c_spn(1)=+0.148816447901042;
    c_spn(2)=-0.132385865767784;
    c_spn(3)=+0.067307604692185;
    c_spn(4)=+0.432666402578175;
    c_spn(5)=0.5-(c_spn(1)+c_spn(2)+c_spn(3)+c_spn(4));
    c_cpl(1)=+0.0502627644003922;
    c_cpl(2)=+0.413514300428344;
    c_cpl(3)=+0.0450798897943977;
    c_cpl(4)=-0.188054853819569;
    c_cpl(5)=+0.541960678450780;
    c_cpl(6)=1-2*(c_cpl(1)+c_cpl(2)+c_cpl(3)+c_cpl(4)+c_cpl(5));

else

    error('spin order not coded')
end

% find the propagators
if nargout>3 || nargin>4

    load_flag=0;
    filename=[param.scratchdir filesep param.job_id filesep 'prop_t' num2str(trt) 'o' num2str(odr) '.mat'];

    % check if propagators are in the cache record
    if strcmp(param.prop_cache,'carry')

        %P=

    elseif strcmp(param.prop_cache,'store') && exist(filename,'file')

        % load the file
        load(filename,'P'); 

        % change the flag to indicate file was just loaded
        load_flag=1;

    elseif nargin<5

        error('interaction Hamiltonian needed for propagator calculation')

    elseif iscell(L1_in)
        
        % default should be the same size as drifts
        P=L1_in;

        for n=1:numel(L1_in)

            dt_trt=dt/trt;
            if odr==0

                P{n}={1}; % don't even use matrices

            elseif odr==1

                P{n}={1,... % no right hand side splitting, don't even use matrices
                    sparse(propagate(param.space,L1_in{n},c_cpl(2)*dt_trt,param.prop_method,param.prop_zeroed))};

            else

                % preallocate cell
                if trt==1

                    % Trotter number 1 has no concatination
                    P{n}=cell(1,length(c_cpl));
                else

                    % Trotter number >1 has concatination of first and last
                    P{n}=cell(1,length(c_cpl)+1);
                end

                % loop over the cpl-constants
                for k=1:length(c_cpl)

                    % ensure a cleanup and sparse propagator
                    P{n}{k}=sparse(propagate(param.space,L1_in{n},c_cpl(k)*dt_trt,param.prop_method,param.prop_zeroed));
                end

                % concatination of first and last if Trotter number >1
                if trt~=1

                    % multiply the first and last propagators
                    P{n}{k+1}=sparse(propagate(param.space,L1_in{n},c_cpl(i_int(1))*dt_trt,param.prop_method,param.prop_zeroed)*...
                        propagate(param.space,L1_in{n},c_cpl(i_int(end))*dt_trt,param.prop_method,param.prop_zeroed));
                end
            end
        end
    else

        dt_trt=dt/trt;
        if odr==0

            P={1}; % don't even use matrices

        elseif odr==1

            P={1,... % no right hand side splitting, don't even use matrices
                sparse(propagate(param.space,L1_in,c_cpl(2)*dt_trt,param.prop_method,param.prop_zeroed))};

        else

            % preallocate cell
            if trt==1

                % Trotter number 1 has no concatination
                P=cell(1,length(c_cpl));
            else

                % Trotter number >1 has concatination of first and last
                P=cell(1,length(c_cpl)+1);
            end

            % loop over the cpl-constants
            for k=1:length(c_cpl)

                % ensure a cleanup and sparse propagator
                P{k}=sparse(propagate(param.space,L1_in,c_cpl(k)*dt_trt,param.prop_method,param.prop_zeroed));
            end

            % concatination of first and last if Trotter number >1
            if trt~=1

                % multiply the first and last propagators
                P{k+1}=sparse(propagate(param.space,L1_in,c_cpl(i_int(1))*dt_trt,param.prop_method,param.prop_zeroed)*...
                    propagate(param.space,L1_in,c_cpl(i_int(end))*dt_trt,param.prop_method,param.prop_zeroed));
            end
        end
    end

    % check if propagators are in the cache record
    if strcmp(param.prop_cache,'store') && ~load_flag
        save(filename,'P','-v7.3');
    end

end

end
