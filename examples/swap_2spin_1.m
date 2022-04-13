clear

% add some paths (relative to this file), so we can find code, you probably
% shouldn't change this if you store all your custom code in this directory
addpath(genpath('../kernel/'))
addpath(genpath('../utilities/'))

% working in a Liouville space (aka adjoint space) with the default
% spherical tensor basis of QOALA functionality

unit  = [+1  0  0  0].';
rho_z = [ 0  0 +1  0].';    rho_z=rho_z/2;

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

% zz product operator
Lzz_12=kron(Jz,unit*rho_z'+rho_z*unit')+kron(unit*rho_z'+rho_z*unit',Jz);

% contruct swap gate from Linden-1999 eq.14
U_swap= (expm(-1i*(+Hz1+Hz2)*3*pi/2)*...
        (expm(-1i*(-Hx1)*pi/2)*...
        (expm(-1i*(+2*Lzz_12)*pi/2)*...
        (expm(-1i*(-Hy1-Hy2)*pi/2)*...
        (expm(-1i*(+2*Lzz_12)*pi/2)*...
        (expm(-1i*(+Hx1+Hx2)*pi/2)*...
        (expm(-1i*(+2*Lzz_12)*pi/2)*...
        (expm(-1i*(+Hy1)*pi/2)))))))));

% clean to swap gate
prop_1=sparse(1e-9*round((1/1e-9)*(U_swap)));

pulse_amplitudes=[1000 1000];

% find the optimised pulse
optpulse=universal_gate_xy('target',prop_1,'interaction',H,'cartops',{ss_ops},...    % required
                        'omega',[0 0],'amplitudes',pulse_amplitudes,...              % required
                        'duration',0.012,'increments',60);                       % required

% calculate the fidelity (ignore the code)
ctrl_sys.pulse_dt=(0.012/60)*ones(1,60);
ctrl_sys.step_method='krylov';
ctrl_sys.space='liouville';
ctrl_sys.ctrl_type='UR';
ctrl_sys.dim=length(H);
drift_sys.drift=H;
disp('... calculating the fidelity %')
[~,fidelity]=waveform_fidelity(ctrl_sys,drift_sys,{Hx1,Hy1,Hx2,Hy2},...
                               2*pi*kron(pulse_amplitudes,[1,1]),optpulse,...
                               speye(size(prop_1)),prop_1);
disp(['fidelity calculated to be ' num2str(100*fidelity,'%.8f') '%'])


% plot the optimal pulse
figure(1)
stairs(0:0.0002:0.012,kron(pulse_amplitudes,[1,1]).*[optpulse; optpulse(end,:)])
xlabel('time (seconds)'); ylabel('pulse amplitude (Hertz)');
legend('x-control 1','y-control 1','x-control 2','y-control 2')
