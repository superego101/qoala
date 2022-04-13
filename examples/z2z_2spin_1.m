clear

% add some paths (relative to this file), so we can find code, you probably
% shouldn't change this if you store all your custom code in this directory
addpath(genpath('../kernel/'))
addpath(genpath('../utilities/'))

% working in a Liouville space (aka adjoint space) with the default
% spherical tensor basis of QOALA functionality
rho_z1=0.5*[0; 0; 0; 0; 0; 0; 0; 0; 1; 0; 0; 0; 0; 0; 0; 0]; % z-magnetisation on spin 1
rho_z2=0.5*[0; 0; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]; % z-magnetisation on spin 2

% this is an interaction Hamiltonian. Here, 2 spins are coupled through,
% what is known as 'the weak coupling approximation' in the language of
% magnetic resonance i.e. only coupled through H_zz terms. The 140 below is
% the coupling strength in Hertz.
H=pi*140* [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;...
            0  0  0  0  0  0  0  0  0 +1  0  0  0  0  0  0;...
            0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;...
            0  0  0  0  0  0  0  0  0  0  0 -1  0  0  0  0;...
            0  0  0  0  0  0 +1  0  0  0  0  0  0  0  0  0;...
            0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;...
            0  0  0  0 +1  0  0  0  0  0  0  0  0  0  0  0;...
            0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;...
            0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;...
            0 +1  0  0  0  0  0  0  0  0  0  0  0  0  0  0;...
            0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;...
            0  0  0 -1  0  0  0  0  0  0  0  0  0  0  0  0;...
            0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0;...
            0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;...
            0  0  0  0  0  0  0  0  0  0  0  0 -1  0  0  0;...
            0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0];

% single spin-1/2 Cartesian rotation operators
I  = [+1  0  0  0 ; 0 +1  0  0 ; 0  0 +1  0 ; 0  0  0 +1];
Jx = [ 0  0  0  0 ; 0  0 +1  0 ; 0 +1  0 +1 ; 0  0 +1  0];  Jx=Jx/sqrt(2);
Jy = [ 0  0  0  0 ; 0  0 -1i 0 ; 0 +1i 0 -1i; 0  0 +1i 0];  Jy=Jy/sqrt(2);
Jz = [ 0  0  0  0 ; 0 +1  0  0 ; 0  0  0  0 ; 0  0  0 -1];

% 2-spin rotation operators of composite system, built from kroneckers of
% 1-spin operators to save space in this example file.
Hx1=sparse(kron(Jx,I)); Hy1=sparse(kron(Jy,I)); Hz1=sparse(kron(Jz,I));
Hx2=sparse(kron(I,Jx)); Hy2=sparse(kron(I,Jy)); Hz2=sparse(kron(I,Jz));

% we need the Cartesian operators for each spin. The arrangement of these
% operators in the cell below is important: the number of rows of the cell
% should be equal to the number of spins; the number of columns depends on
% how we control the spins. Below is the example of each of the two spins
% having their individual controls
ss_ops={Hx1,Hy1,Hz1,[],[],[];
        [],[],[],Hx2,Hy2,Hz2};

pulse_amplitudes=[1000 1000];

% find the optimised pulse
optpulse=state2state_xy('initial',rho_z1,'target',rho_z2,'interaction',H,'cartops',{ss_ops},...    % required
                        'omega',[0 0],'amplitudes',pulse_amplitudes,...              % required
                        'duration',0.01,'increments',50);                       % required

% calculate the fidelity (ignore the code)
ctrl_sys.pulse_dt=(0.01/50)*ones(1,50);
ctrl_sys.step_method='krylov';
ctrl_sys.space='liouville';
ctrl_sys.ctrl_type='PP';
drift_sys.drift=H;
disp('... calculating the fidelity %')
[~,fidelity]=waveform_fidelity(ctrl_sys,drift_sys,{Hx1,Hy1,Hx2,Hy2},...
                               2*pi*kron(pulse_amplitudes,[1,1]),optpulse,...
                               rho_z1/norm(rho_z1),rho_z2/norm(rho_z2));
disp(['fidelity calculated to be ' num2str(100*fidelity,'%.8f') '%'])


% plot the optimal pulse
figure(1)
stairs(0:0.0002:0.01,kron(pulse_amplitudes,[1,1]).*[optpulse; optpulse(end,:)])
xlabel('time (seconds)'); ylabel('pulse amplitude (Hertz)');
legend('x-control 1','y-control 1','x-control 2','y-control 2')
