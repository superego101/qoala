clear

% add some paths (relative to this file), so we can find code, you probably
% shouldn't change this if you store all your custom code in this directory
addpath(genpath('../kernel/'))
addpath(genpath('../utilities/'))

% working in a Liouville space (aka adjoint space) with the default
% spherical tensor basis of QOALA functionality
unit  = [+1  0  0  0].';
rho_z = [ 0  0 +1  0].';    rho_z=rho_z/2;
rho_z1=kron(kron(rho_z,unit),unit);             % z-magnetisation on spin 1
rho_z3=kron(unit,kron(unit,rho_z));             % z-magnetisation on spin 3

% single spin-1/2 Cartesian rotation operators
I  = [+1  0  0  0 ; 0 +1  0  0 ; 0  0 +1  0 ; 0  0  0 +1];
Jx = [ 0  0  0  0 ; 0  0 +1  0 ; 0 +1  0 +1 ; 0  0 +1  0];  Jx=Jx/sqrt(2);
Jy = [ 0  0  0  0 ; 0  0 -1i 0 ; 0 +1i 0 -1i; 0  0 +1i 0];  Jy=Jy/sqrt(2);
Jz = [ 0  0  0  0 ; 0 +1  0  0 ; 0  0  0  0 ; 0  0  0 -1];

% this is an interaction Hamiltonian. Here, 3 spins are coupled in a linear 
% chain through, what is known as 'the weak coupling approximation' in the 
% language of magnetic resonance i.e. only coupled through H_zz terms. The 
% 140 and -160 below are the coupling strengths in Hertz.
Lzz_12=kron(kron(Jz,unit*rho_z'+rho_z*unit'),I)+kron(kron(unit*rho_z'+rho_z*unit',Jz),I);
Lzz_23=kron(I,kron(Jz,unit*rho_z'+rho_z*unit'))+kron(I,kron(unit*rho_z'+rho_z*unit',Jz));
H=2*pi*(140)*Lzz_12 + 2*pi*(-160)*Lzz_23;

% 3-spin rotation operators of composite system, built from kroneckers of
% 1-spin operators to save space in this example file.
Hx1=sparse(kron(kron(Jx,I),I)); Hy1=sparse(kron(kron(Jy,I),I)); Hz1=sparse(kron(kron(Jz,I),I));
Hx2=sparse(kron(kron(I,Jx),I)); Hy2=sparse(kron(kron(I,Jy),I)); Hz2=sparse(kron(kron(I,Jz),I));
Hx3=sparse(kron(I,kron(I,Jx))); Hy3=sparse(kron(I,kron(I,Jy))); Hz3=sparse(kron(I,kron(I,Jz)));

% we need the Cartesian operators for each spin. The arrangement of these
% operators in the cell below is important: the number of rows of the cell
% should be equal to the number of spins; the number of columns depends on
% how we control the spins. Below is the example of each of the two spins
% having their individual controls
ss_ops={Hx1,Hy1,Hz1,[],[],[],[],[],[];
        [],[],[],Hx2,Hy2,Hz2,[],[],[];
        [],[],[],[],[],[],Hx3,Hy3,Hz3};

pulse_amplitudes=[1000 1000 1000];

% find the optimised pulse
optpulse=state2state_xy('initial',rho_z1,'target',rho_z3,'interaction',H,'cartops',{ss_ops},...    % required
                        'omega',[0 0 0],'amplitudes',pulse_amplitudes,...              % required
                        'duration',0.022,'increments',110);                       % required

% calculate the fidelity (ignore the code)
ctrl_sys.pulse_dt=(0.022/110)*ones(1,110);
ctrl_sys.step_method='krylov';
ctrl_sys.space='liouville';
ctrl_sys.ctrl_type='PP';
drift_sys.drift=H;
disp('... calculating the fidelity %')
[~,fidelity]=waveform_fidelity(ctrl_sys,drift_sys,{Hx1,Hy1,Hx2,Hy2,Hx3,Hy3},...
                               2*pi*kron(pulse_amplitudes,[1,1]),optpulse,...
                               rho_z1/norm(rho_z1),rho_z3/norm(rho_z3));
disp(['fidelity calculated to be ' num2str(100*fidelity,'%.8f') '%'])


% plot the optimal pulse
figure(1)
stairs(0:(0.022/110):0.022,kron(pulse_amplitudes,[1,1]).*[optpulse; optpulse(end,:)])
xlabel('time (seconds)'); xlim([0 0.022]); ylabel('pulse amplitude (Hertz)');
legend('x-control 1','y-control 1','x-control 2','y-control 2','x-control 3','y-control 3')
