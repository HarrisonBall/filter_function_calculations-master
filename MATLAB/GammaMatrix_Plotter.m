function [PLOT_GAMMA_MATRIX_CONTROL_SEQUENCE] = GammaMatrix_Plotter(GammaMatrix)

%--------------------------------------------------------------------------
%                                   TEST 
%--------------------------------------------------------------------------
% clear all
% clf('reset')
% pitime = 1;
% theta = pi;
% GammaMatrix = GammaMatrix_MultiAxis(pitime);
% GammaMatrix_Plotter(GammaMatrix)

%__________________________________________________________________________

%                               REFERENCES 
%__________________________________________________________________________

%   [1] H. Ball and M. J. Biercuk, "Walsh-synthesized noise-filtering 
%       quantum logic", arXiv:1410.1624, EPJ Quantum Technology 2, 1 (2015).
%   [2] Todd J. Green, Jarrah Sastrawan, Hermann Uys, and M.J. Biercuk, 
%       "Arbitrary quantum control of qubits in the presence of universal 
%       noise."  New Journal of Physics 15, 095004 (2013).
%__________________________________________________________________________

%                            FUNCTION DESCRIPTION
%__________________________________________________________________________

%This function generates a plot of the various control fields describing 
%a time-dependent control Hamiltonian for a GENERAL CONTROL setting (see
%bellow). The control space of said Hamiltonian must be encoded as a 
%standard GammaMatrix. 
%__________________________________________________________________________

%                     STANDARD GAMMA-MATRIX NOTATION
%__________________________________________________________________________

%The standardized GammaMatrix parameterization of a control sequence is set
%out in detail in the supporting mathematica notebook 
%[H.Ball 2016] filter_function_calculator (MASTER).nb. 
%This parameterization is consistent with (though generalizes) the notation 
%introudced in Eq. 11 of Ref. [1]. This framework generically describes any 
%arbitrary N-segment unitary control sequence arising from a piecewise-constant, 
%time-dependent control Hamiltonian (see Eq. 1 of [1] & Eq. 5 of [2]). 

%Such a control Hamiltonian generally results a sequence of N arbitrary, 
%finite-duration rotations of the Bloch vector, corresponding to a sequence 
%of N (arbitary) unitary control operations on a single-qubit. This is 
%referred to as GENERAL CONTROL. Each pulse segment in the sequence is 
%parameterized by a Rabi rate, pulse duration, angle of rotation, and 
%direction of rotation. For GENERAL CONTROL, the direction of rotation 
%is specifed by a unit vector, nhat, pointing in any direction in the Bloch
%sphere:

%           GENERAL CONTROL:           nhat = (nx, ny, nz).

%That is, each unitary may enact a distinct rotation of the Bloch vector
%about any axis in R^3. A common special case is 

%           EQUATORIAL CONTROL:        nhat = (cos(phi), sin(phi), 0),

%That is, the direction of control is restricted to the equatorial plane 
%of the Bloch sphere (nz = 0), so the direction of the rotation axis may be
%parameterized by a "phase" angle, phi. A further special case is 

%           SINGLE-AXIS (X) CONTROL:   nhat = (1, 0, 0)

%corresponding to the choice phi = 0. 

%--------------------------------------------------------------------------
%                            CRITICAL FORM! 
%--------------------------------------------------------------------------

%The form of the standard GammaMatrix is critical! The required form of the 
%output of this function (and any other GammaMatrix_[text].m file) is a 
%6-column matrix of the form:  

%              (Omegal)  (taul) (theta)    ..... (nhat) .....               

%                 c1      c2      c3       c4      c5     c6            

%             __________________________________________________

%     r1      | Omega1 | tau1 | theta1 |  n1x  |  n1y  |  n1z  |
%             __________________________________________________

%     r2      | Omega2 | tau2 | theta2 |  n2x  |  n2y  |  n2z  |  
%             __________________________________________________

%     r3      | Omega3 | tau3 | theta3 |  n3x  |  n3y  |  n3z  |  
%             __________________________________________________

%                               .           .
%                               .           .
%             __________________________________________________

%     rl      | Omegal | taul | thetal |  nlx  |  nly  |  nlz  |  
%             __________________________________________________

%                               .           .
%                               .           .
%             __________________________________________________

%     rN      | OmegaN | tauN | thetaN |  nNx  |  nNy  |  nNz  |
%             __________________________________________________

%where we define


%     N        ...       total number of pulses in the composite sequence.

%     l        ...       index labelling pulses in sequence. This index
%                        runs from 1 to N.

%     Omegal   ...       Rabi rate during the lth pulse.

%     taul     ...       duration of lth pulse (in seconds).

%     thetal   ...       angle of rotation swept out during the lth pulse 
%                        (in radians). Note: thetal = Omegal*taul. This is
%                        the angle through which the Bloch vector rotates
%                        under the lth control pulse segment. 

%     nhat     ...       unit vector (nlx,nly,nlz) defining rotation axis 
%                        during the lth pulse.

%     nlx      ...       x-component of rotation axis during lth pulse.
%     nly      ...       y-component of rotation axis during lth pulse.
%     nlz      ...       z-component of rotation axis during lth pulse.
  
%__________________________________________________________________________

%                      UNPACK GammaMatrix DATA
%__________________________________________________________________________

micro = 1E-6;
G = GammaMatrix;

% G(:,1) = Omega_col;
% G(:,2) = tau_col;
% G(:,3) = theta_col;
% G(:,4) = nx_col;
% G(:,5) = ny_col;
% G(:,6) = nz_col;

%Total duration of control sequence (dimensionalised, in seconds)
tau_tot = sum(G(:,2));

%Rotation Axis Components

nx_col = G(:,4);
ny_col = G(:,5);
nz_col = G(:,6);

%Pulse Durations (dimensionalised, in seconds)
PulseDurations_Dimensionalised = G(:,2);

%Fractional Pulse Durations (dimensionalised pulse durations/tau_tot)
%PulseDurations_Normalized = G(:,2)/tau_tot;

%Rabi Rates (dimensionalised)
RabiRates_Dimensionalised = G(:,1);

%Bloch Rotation Angles 
BlochRotationAngles = G(:,3);

%Total Number of Pulse Segments in Sequence
N = length(G(:,1));

%Pulse Segment Start/End Times (dimensionalised, in seconds)
PulseBoundaries_Dimensionalised = zeros(N+1,1);

for i=1:N
    PulseBoundaries_Dimensionalised(i+1) = sum(PulseDurations_Dimensionalised(1:i));
end

%__________________________________________________________________________

%                                 PLOT
%__________________________________________________________________________


%net gate angle (NGA)
NGA = mod(sum(BlochRotationAngles),2*pi);


%--------------------------------------------------------------------------
%                                SETUP
%--------------------------------------------------------------------------

%Rotation Axis Components (for stair plotter) 
nx_col_StairPlotter = nx_col;
nx_col_StairPlotter(N+1) = nx_col(end);

ny_col_StairPlotter = ny_col;
ny_col_StairPlotter(N+1) = ny_col(end);

nz_col_StairPlotter = nz_col;
nz_col_StairPlotter(N+1) = nz_col(end);

%Rotation Axis Components (for area plotter) 

nx_col_area_plot = zeros(1,2*N);
nx_col_area_plot(1:2:2*N) = nx_col;
nx_col_area_plot(2:2:2*N) = nx_col;

ny_col_area_plot = zeros(1,2*N);
ny_col_area_plot(1:2:2*N) = ny_col;
ny_col_area_plot(2:2:2*N) = ny_col;

nz_col_area_plot = zeros(1,2*N);
nz_col_area_plot(1:2:2*N) = nz_col;
nz_col_area_plot(2:2:2*N) = nz_col;

%Rabi Rates (for stair plotter) 
RabiRates_StairPlotter = RabiRates_Dimensionalised;
RabiRates_StairPlotter(N+1) = RabiRates_Dimensionalised(end);

%Bloch Rotation Angles (for stair plotter) 
BlochRotationAngles_StairPlotter = BlochRotationAngles;
BlochRotationAngles_StairPlotter(N+1) = BlochRotationAngles(end);

%Normalized Pulse Segment Stard/End Times 
PulseBoundaries_NonDim = PulseBoundaries_Dimensionalised/tau_tot;

%Normalized Pulse Segment Stard/End times (for area plotter) 
PulseBoundaries_NonDim_area_plot= zeros(1,2*N);
PulseBoundaries_NonDim_area_plot(1:2:2*N) = PulseBoundaries_NonDim(1:N);
PulseBoundaries_NonDim_area_plot(2:2:2*N) = PulseBoundaries_NonDim(2:N+1);


%--------------------------------------------------------------------------
%                    SUBPLOT #1: Rotation Angles 
%--------------------------------------------------------------------------
figure(100)
clf('reset')
subplot(5,1,1)

stairs(PulseBoundaries_NonDim,BlochRotationAngles_StairPlotter/pi);

titlestring {1} = ['Gate Angle = ', num2str(NGA/pi), '\pi'];
titlestring {2} = ['Sequence Duration \tau = ', num2str(tau_tot/micro),' \mus'];
title(titlestring)


hold on
stem(PulseBoundaries_NonDim,BlochRotationAngles_StairPlotter/pi,':r','Marker','none');

ylabel('Rotation Angle \theta (\pi)');

%Axis scaling
if min(BlochRotationAngles)*max(BlochRotationAngles)~=0;
    BRA_clearance = 0.2*range(BlochRotationAngles_StairPlotter/pi);
    BRA_ymin = min([0 min(BlochRotationAngles_StairPlotter/pi)-BRA_clearance;]);    
    BRA_ymax = max([max(BlochRotationAngles_StairPlotter/pi)+BRA_clearance,...
        1.2*max(BlochRotationAngles_StairPlotter/pi)]);
   
    axis([0,max(PulseBoundaries_NonDim),BRA_ymin,BRA_ymax]);
end

%--------------------------------------------------------------------------
%                    SUBPLOT #2: Rabi Rates
%--------------------------------------------------------------------------
subplot(5,1,2)

stairs(PulseBoundaries_NonDim,RabiRates_StairPlotter/pi);
hold on
stem(PulseBoundaries_NonDim,RabiRates_StairPlotter/pi,':r','Marker','none');

ylabel('Rabi Rate \Omega (\pi)');

%Axis scaling
if min(RabiRates_StairPlotter/pi)*max(RabiRates_StairPlotter/pi)~=0;
       Rabi_clearance = 0.2*range(RabiRates_StairPlotter/pi);
    Rabi_ymin = min([0 min(RabiRates_StairPlotter/pi)-Rabi_clearance;]);    
    Rabi_ymax = max([max(RabiRates_StairPlotter/pi)+Rabi_clearance,...
        1.2*max(RabiRates_StairPlotter/pi)]);
   
    axis([0,max(PulseBoundaries_NonDim),Rabi_ymin,Rabi_ymax]);
end

%--------------------------------------------------------------------------
%              SUBPLOT #3: nx (rotation axis component)   
%--------------------------------------------------------------------------

subplot(5,1,3)

area(PulseBoundaries_NonDim_area_plot,nx_col_area_plot,'FaceColor',[0,0,0.6],'EdgeColor','none');
hold on
stairs(PulseBoundaries_NonDim,nx_col_StairPlotter,':white')
hold on
stem(PulseBoundaries_NonDim,nx_col_StairPlotter,':white','Marker','none');
axis([0,1,-1.1,1.1])
ylabel('n_x');

%--------------------------------------------------------------------------
%              SUBPLOT #3: ny (rotation axis component)   
%--------------------------------------------------------------------------

subplot(5,1,4)

area(PulseBoundaries_NonDim_area_plot,ny_col_area_plot,'FaceColor',[0,0,0.6],'EdgeColor','none');
hold on
stairs(PulseBoundaries_NonDim,ny_col_StairPlotter,':white')
hold on
stem(PulseBoundaries_NonDim,ny_col_StairPlotter,':white','Marker','none');
axis([0,1,-1.1,1.1])
ylabel('n_y');

%--------------------------------------------------------------------------
%              SUBPLOT #3: nz (rotation axis component)   
%--------------------------------------------------------------------------

subplot(5,1,5)

area(PulseBoundaries_NonDim_area_plot,nz_col_area_plot,'FaceColor',[0,0,0.6],'EdgeColor','none');
hold on
stairs(PulseBoundaries_NonDim,nz_col_StairPlotter,':white')
hold on
stem(PulseBoundaries_NonDim,nz_col_StairPlotter,':white','Marker','none');
axis([0,1,-1.1,1.1])
xlabel('time (\tau)');
ylabel('n_z');

end



