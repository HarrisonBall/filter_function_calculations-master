

function  LOCAL_AMPLITUDE_CONTROL_VECTOR_EQUATORIAL_CONTROL_GAMMA_MATRIX = Local_AmplitudeControlVector_EquatorialControl_GammaMatrix(GammaMatrix,l)

%--------------------------------------------------------------------------
%                                   TEST 
%--------------------------------------------------------------------------
% 
% clear all
% format long
% pitime = 1/4;
% w_col = 0.1;
% GammaMatrix = GammaMatrix_DCG(pitime);
% l=1;
% Local_AmplitudeControlVector_EquatorialControl_GammaMatrix(GammaMatrix,l)

%__________________________________________________________________________

%                              REFERENCES 
%__________________________________________________________________________

%   [1] H. Ball and M. J. Biercuk, "Walsh-synthesized noise-filtering 
%       quantum logic", arXiv:1410.1624, EPJ Quantum Technology 2, 1 (2015).
%   [2] Todd J. Green, Jarrah Sastrawan, Hermann Uys, and M.J. Biercuk, 
%       "Arbitrary quantum control of qubits in the presence of universal 
%       noise."  New Journal of Physics 15, 095004 (2013).

%__________________________________________________________________________

%                         FUNCTION DESCRIPTION
%__________________________________________________________________________

%This function outputs the (implicitely w-domain) "local amplitude control 
%vector" (LACV) for the "lth" pulse, associated (multuplicative) amplitude 
%damping noise coaxial with the direction of EQUATORIAL control. That is, 
%if lth control pulse enacts a rotation about the equatorial axis 
%sig_phi = cos(phi)sig_x + sin(phi)sig_y then, during this time interval, 
%the noise is also applied along the sig_phi axis. This model is
%appropriate for considering muluplicative noise in the amplitude of a 
%driving field, manifesting in time-dependent fluctuations proportional to 
%the Rabi rate. in the Rabi rate (see Eq. 3 of Ref. [1]). The form of the 
%LACV can be seen in Eq. 16 of Ref. [1]. Namely, the "l-1th" Control 
%History Matrix left multiplied by the "lth" Segment Projection Vector. 

%The control Hamiltonian is assumed to satisfy the EQUATORIAL control form,
%a special case of GENERAL control, described below. 

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

%                            FUNCTION INPUTS 
%__________________________________________________________________________

%The required inputs for this m-file are as follows: 

%     GammaMatrix        ...       standard GammaMatrix parameterizing the
%                                  piecewise-constant control sequence. 
%     l                  ...       index labelling the "lth" pulse in the
%                                  sequence. This index runs from 
%                                  1 to N, where N = length(GammaMatrix) is 
%                                  the total number of pulse segments in
%                                  the control sequence. 

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

%                            FUNCTION OUTPUT 
%__________________________________________________________________________

% 1x3 vector, consisting of the wieghted sum of the first two rows
% of the l-1th History Matrix. The wieghts are given by the elements of
% the lth Projection Vector. 

%__________________________________________________________________________

%                          UNPACK GammaMatrix DATA
%__________________________________________________________________________

G = GammaMatrix;

% G(:,1) = Omega_col;
% G(:,2) = tau_col;
% G(:,3) = theta_col;
% G(:,4) = nx_col;
% G(:,5) = ny_col;
% G(:,6) = nz_col;

%Rotation Axis Components of lth Pulse Segment
nx_l = G(l,4);
ny_l = G(l,5);
%nz_L = G(l,6);     %set to zero since we ssume EQUATORIAL control

%Rabi Rate of lth Pulse Segment (dimensionalised)
Omega_l = G(l,1);

%__________________________________________________________________________

%                      DEFINE PROJECTION VECTOR
%__________________________________________________________________________

ProjectionVector_l = (Omega_l/2)*[nx_l ny_l 0];

%__________________________________________________________________________

%                             OUTPUT
%__________________________________________________________________________

LOCAL_AMPLITUDE_CONTROL_VECTOR_EQUATORIAL_CONTROL_GAMMA_MATRIX=ProjectionVector_l*HistoryMatrix_GammaMatrix(GammaMatrix,l-1);

end


