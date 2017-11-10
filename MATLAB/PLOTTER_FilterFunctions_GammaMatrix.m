%__________________________________________________________________________

%                            INITIALIZE SCRIPT 
%__________________________________________________________________________

clear all
clf('reset')

%__________________________________________________________________________

%                              REFERENCES 
%__________________________________________________________________________

%   [1] H. Ball and M. J. Biercuk, "Walsh-synthesized noise-filtering 
%       quantum logic", arXiv:1410.1624, EPJ Quantum Technology 2, 1 (2015).
%   [2] Todd J. Green, Jarrah Sastrawan, Hermann Uys, and M.J. Biercuk, 
%       "Arbitrary quantum control of qubits in the presence of universal 
%       noise."  New Journal of Physics 15, 095004 (2013).

%__________________________________________________________________________

%                         SCRIPT DESCRIPTION
%__________________________________________________________________________

%This script plots the w-domain filter function associated with the 
%noise in the dephasing (sigma_z) quadrature, as set out in Eq. 13 of 
%Ref. [1]. The control Hamiltonian is assumed to satisfy the GENERAL 
%control form. 

%As set out in the supporting mathematica notebooks, GENERAL control
%assumes a control Hamiltonian implementing a sequence of N arbitrary 
%finite-duration rotations of the Bloch vector, corresponding to a sequence 
%of N (arbitary) unitary control operations on a single-qubit. Each pulse
%segment in the sequence is parameterized by a Rabi rate, pulse duration, 
%angle of rotation, and direction of rotation. The direction of rotation 
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

%__________________________________________________________________________

%                            SCRIPT INPUTS 
%__________________________________________________________________________

%The required inputs for this m-file are as follows: 

%     GammaMatrix        ...       standard GammaMatrix parameterizing the
%                                  piecewise-constant control sequence. 

%     w_col              ...       column vector of ANGULAR frequencies.

%--------------------------------------------------------------------------
%                            CRITICAL FORM! 
%--------------------------------------------------------------------------

%The form of the standard GammaMatrix is critical! This form follows 
%the same arrangement presented in the supporting mathematica notebooks:

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

%                            	WARNING!
%__________________________________________________________________________

%Care should be taken that the right m.file is used to compute 
%the filter function (or other objects), depending on which of the above 
%regimes the in put GammaMatrix falls into. For instance, a filter function 
%computation may be sped up significantly, if one produces a computation 
%specialized to single-axis control. However, if one uses this single-axis 
%specialized package to compute the filter function for an input GammaMatrix 
%parameterizing a control Hamiltonian of the GENERAL control variety, the 
%result will be nonsenseical!

%__________________________________________________________________________

%                DEFINE TIME/FREQUENCY DOMAIN PARAMETERS   
%__________________________________________________________________________

pitime = 1000E-6;

w_IndexStart_Normalized = -6;
w_IndexEnd_Normalized = 1;

w_IndexStart_InvPitimeScaled = w_IndexStart_Normalized - log10(pitime);
w_IndexEnd_InvPitimeScaled = w_IndexEnd_Normalized - log10(pitime);

w_col = logspace(w_IndexStart_InvPitimeScaled,w_IndexEnd_InvPitimeScaled,2000)'*2*pi;

%__________________________________________________________________________

%                          GAMMA MATRIX LIBRARY
%__________________________________________________________________________

%G = GammaMatrix_PRIM(pi,pitime);
G = GammaMatrix_DCG(pitime);
%G = GammaMatrix_BB1(pi/2,pitime);

%__________________________________________________________________________

%                 COMPUTE FILTER FUNCTIONS (UNMODIFIED)
%__________________________________________________________________________

%Quadrature: dephasing (sigma_z noise)
FF_z = FF_GammaMatrix_GeneralControl_Dephasing(G,w_col);
%Quadrature: stationary x-relaxiation (sigma_x noise)
FF_x = FF_GammaMatrix_GeneralControl_RelaxationX(G,w_col);
%Quadrature: coaxial amplitude noise (sigma_phi noise)
FF_O = FF_GammaMatrix_EquatorialControl_MultuplicativeAmplitude(G,w_col);

%__________________________________________________________________________

%                   COMPUTE FILTER FUNCTIONS (MODIFIED)
%__________________________________________________________________________

%Quadrature: dephasing (sigma_z noise)
FF_z_modified = FF_z./(w_col).^2;
%Quadrature: stationary x-relaxiation (sigma_x noise)
FF_x_modified = FF_x./(w_col).^2;
%Quadrature: coaxial amplitude noise (sigma_phi noise)
FF_O_modified = FF_O./(w_col).^2;

%__________________________________________________________________________

%                      UNPACK GammaMatrix DATA
%__________________________________________________________________________

micro = 1E-6;

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
PulseDurations_Normalized = G(:,2)/tau_tot;

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

%                         PLOT FILTER FUNCTIONS
%__________________________________________________________________________

figure(30)
clf('reset')
%--------------------------------------------------------------------------
%                     ANGULAR FREQUENCY DIMENSIONING 
%--------------------------------------------------------------------------

%Define scaling unit for graphing unmodified FFs
w_col_scaling_unit_unmodified = 2*pi/pitime;
%Define scaling unit for graphing modified FFs
w_col_scaling_unit_modified = 2*pi/pitime;

%Rescale w_col in units of the above scaling unit for graphing unmodified FFs
w_col_scaled_unmodified = w_col/w_col_scaling_unit_unmodified;
%Rescale w_col in units of the above scaling unit for graphing modified FFs
w_col_scaled_modified = w_col/w_col_scaling_unit_modified;

%--------------------------------------------------------------------------
%              DEFINE ORDER OF ERROR SUPPRESSION EXPONENTS 
%--------------------------------------------------------------------------

alpha2 = 2;
alpha4 = 4;
alpha6 = 6;

%--------------------------------------------------------------------------
%                              SUBPLOT #1: 
%           Plot unmodified FFs, & order of error suppression 
%--------------------------------------------------------------------------

subplot(2,1,1)

loglog(w_col_scaled_unmodified,FF_z,...
       w_col_scaled_unmodified,FF_O,...    
       w_col_scaled_unmodified,FF_x)       %CHANGE THIS
 
 xlabel('\omega   (2\pi/\tau_{\pi})','FontSize',14);
 ylabel('F(\omega)','FontSize',14);
 
 title('Unmodified Filter Functions','FontSize',14);
 
hold on

loglog(w_col_scaled_unmodified,w_col_scaled_unmodified.^alpha2,'black--',...
   w_col_scaled_unmodified,w_col_scaled_unmodified.^alpha4,'b --',...
   w_col_scaled_unmodified,w_col_scaled_unmodified.^alpha6,'r--')

legend30211 = legend('F_z(\omega)',...
   'F_\Omega(\omega)',...
   'F_x(\omega)',...
   '(p - 1) = 0',...
   '(p - 1) = 1',...
   '(p - 1) = 2',...
   'Location','NorthEastOutside');
set(legend30211,'FontSize',14);

axis([w_col_scaled_unmodified(1) w_col_scaled_unmodified(end) 10^(-35) 10^(10)]) 

%--------------------------------------------------------------------------
%                              SUBPLOT #2: 
%            Plot modified FFs, & order of error suppression 
%--------------------------------------------------------------------------

subplot(2,1,2)

loglog(w_col_scaled_modified,FF_z_modified,...
           w_col_scaled_modified,FF_O_modified,...   
           w_col_scaled_modified,FF_x_modified)      %CHANGE THIS
 
 xlabel('\omega   (2\pi/\tau_{\pi})','FontSize',14);
 ylabel('F(\omega)/\omega^2','FontSize',14);

 title('Modified Filter Functions','FontSize',14);

hold on

loglog(w_col_scaled_modified,w_col_scaled_modified.^(alpha2-2),'black--',...
   w_col_scaled_modified,w_col_scaled_modified.^(alpha4-2),'b --',...
   w_col_scaled_modified,w_col_scaled_modified.^(alpha6-2),'r--')

legend30212 = legend('F_z(\omega)/\omega^2',...
   'F_\Omega(\omega)/\omega^2',...
   'F_x(\omega)/\omega^2',...
   '(p - 1) = 0',...
   '(p - 1) = 1',...
   '(p - 1) = 2',...
   'Location','NorthEastOutside');
set(legend30212,'FontSize',14);

axis([w_col_scaled_modified(1) w_col_scaled_modified(end) 10^(-35) 10^(10)]) 


%--------------------------------------------------------------------------
%                         PLOT control sequence
%--------------------------------------------------------------------------

GammaMatrix_Plotter(G);


