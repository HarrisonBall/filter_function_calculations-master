
function [FF_GAMMA_MATRIX_EQUATORIAL_MULTUPLICATIVE_AMPLITUDE] = FF_GammaMatrix_EquatorialControl_MultuplicativeAmplitude(GammaMatrix,w_col)

%--------------------------------------------------------------------------
%                                   TEST 
%--------------------------------------------------------------------------

% clear all
% format long
% pitime = 1/4;
% w_col = 0.1;
% GammaMatrix = GammaMatrix_DCG(pitime);
% FF_GammaMatrix_EquatorialControl_MultuplicativeAmplitude(GammaMatrix,w_col)

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

%This function outputs the w-domain filter function associated with 
%noise in the (multuplicative) amplitude damping qudrature. That is, noise
%coaxial with the direction of EQUATORIAL control. This model is
%appropriate for considering muluplicative noise in the amplitude of a 
%driving field, manifesting in time-dependent fluctuations proportional to 
%the Rabi rate. in the Rabi rate (see Eq. 3 of Ref. [1]). The form of 
%filter function is set out in Eq. 14 of Ref. [1].

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

%                         COMPUTATIONAL OBJECTS
%__________________________________________________________________________

%Imaginary number
ii=sqrt(-1);

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

%Total number of pulse segments in control sequence
N = length(G(:,1));

%Pulse Durations (dimensionalised, in seconds) 
tau_col = G(:,2);

%__________________________________________________________________________

%           OBTAIN PULSE SEGMENT START/END TIMES: [t_{l-1}, t_{l}]
%__________________________________________________________________________

time = zeros(N+1,1);

for i=1:N
    time(i+1) = sum(tau_col(1:i));
end

%__________________________________________________________________________        

%              NUMERICALLY COMPUTE DEPHASING FILTER FUNCTION  
%__________________________________________________________________________
                                
%--------------------------------------------------------------------------
%                   INITIALIZE SIMULATION QUANTITIES 
%--------------------------------------------------------------------------

R_Omega = zeros(length(w_col),3);

%--------------------------------------------------------------------------
%               COMPUTE TOTAL DEPHASING CONTROL VECTOR  
%                       see Eq. 127 of Ref. [1] 
%--------------------------------------------------------------------------

for l=1:N;                          %Index the n time intervals by l      
    
    EXPtl = exp(ii*time(l+1)*w_col); %NOTE: time(i+1) = t_{(i-1)+1}
    EXPtlm1 = exp(ii*time(l)*w_col); %NOTE: time(i) = t_{i-1}
    
    EXPfactor = EXPtlm1-EXPtl;
    
    A = Local_AmplitudeControlVector_EquatorialControl_GammaMatrix(GammaMatrix,l);
    
    %A = [RabiRate(l) 0 0];
    %Simplification by K.S. True for single-axis amplitude-modulated
    %control. Not general result. 
            
    Update = EXPfactor*A;
    
    R_Omega = R_Omega+Update;
            
end

%--------------------------------------------------------------------------
%                       DEPHASING FILTER FUNCTION 
%--------------------------------------------------------------------------

Filt_Omega = sum(R_Omega.*conj(R_Omega),2);

%____________________________________________________________________
        
%                         FUNCTION OUTPUT
%____________________________________________________________________

FF_GAMMA_MATRIX_EQUATORIAL_MULTUPLICATIVE_AMPLITUDE = Filt_Omega;
%output: column vector. 

end





