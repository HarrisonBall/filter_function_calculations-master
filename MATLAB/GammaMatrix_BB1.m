function [GAMMA_BB1] = GammaMatrix_BB1(theta,pitime)

%--------------------------------------------------------------------------
%                                   TEST 
%--------------------------------------------------------------------------
% clear all
% pitime = 1/4;
% theta = pi/2;
% GammaMatrix_BB1(pitime,theta)

%__________________________________________________________________________

%                               REFERENCES 
%__________________________________________________________________________

%   [1] H. Ball and M. J. Biercuk, "Walsh-synthesized noise-filtering 
%       quantum logic", arXiv:1410.1624, EPJ Quantum Technology 2, 1 (2015).
%   [2] Todd J. Green, Jarrah Sastrawan, Hermann Uys, and M.J. Biercuk, 
%       "Arbitrary quantum control of qubits in the presence of universal 
%       noise."  New Journal of Physics 15, 095004 (2013).
%__________________________________________________________________________

%                         FUNCTION DESCRIPTION
%__________________________________________________________________________

%This m-file outputs the standard GammaMatrix-parametrization of a BB1 
%composite pulse sequence (Eq. 44 of Ref. [1]). In this implementation, we
%assume constant Rabi rate, and variable pulse segment durations. 
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

%                                OUTPUT
%__________________________________________________________________________

%Rabi Rates
Omega_col = (pi/pitime)*ones(4,1);

%Rotation Angles
theta_col = [theta, pi, 2*pi,pi]';

%Rotation Axes
phi_1 = acos(-theta/(4*pi));

nhat_col = [
    1,0,0;
    cos(phi_1), sin(phi_1),0;
    cos(3*phi_1), sin(3*phi_1),0;
    cos(phi_1), sin(phi_1),0;
    ];

tau_col = theta_col./Omega_col;
  
GAMMA_BB1  = cat(2,Omega_col,tau_col,theta_col,nhat_col);
 
end





