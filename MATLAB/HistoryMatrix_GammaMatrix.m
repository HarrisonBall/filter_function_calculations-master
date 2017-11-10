

function HISTORY_MATRIX_GAMMA_MATRIX = HistoryMatrix_GammaMatrix(GammaMatrix,k)

%--------------------------------------------------------------------------
%                                   TEST 
% --------------------------------------------------------------------------
% clear all
% clf('reset')
% pitime = 1;
% GammaMatrix = GammaMatrix_PRIM(pi,pitime);
% tau = pitime;
% k=0;
% HistoryMatrix_GammaMatrix(GammaMatrix,k)

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

%This function outputs the Control History Matrix (see Eq. 80 in 
%Ref. [1]) associated with the "kth" pulse segment of a piecewise-constant 
%control pulse sequence. The control Hamiltonian is assumed to satisfy the
%GENERAL control form. 

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

%                            FUNCTION INPUTS 
%__________________________________________________________________________

%The required inputs for this m-file are as follows: 

%     k                 ...        index labelling the number of completed 
%                                  pulses. That is, HistoryMatrix(k) 
%                                  captures the result of all pulse segments
%                                  1,2,...,k (inclusive of k). This index
%                                  runs from 1 through to N-1, where 
%                                  N = length(GammaMatrix) is the total 
%                                  number of pulse segments in the control 
%                                  sequence. 

%     GammaMatrix        ...       standard GammaMatrix parameterizing the
%                                  piecewise-constant control sequence. 

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

% 3x3 matrix, whose ijth element is given by 
% 0.5*Trace(Q^{(k)}^dagger*sigma_j*Q^{(k)}sigma_j), where Q is the
% Cumulative Operator defined in Eq. 9 of Ref. [1].

%__________________________________________________________________________

%                         COMPUTATIONAL OBJECTS
%__________________________________________________________________________

%Imaginary number
ii=sqrt(-1);

%Define Pauli Matrices 
I=eye(2,2);
sigx=[0 1 ; 1 0];
sigy=[0 -ii; ii 0];
sigz=[1 0; 0 -1];

%Initialize 
Q=I;

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

%tau_tot = sum(G(:,2));

%Rotation Axis Components
nx_col = G(:,4);
ny_col = G(:,5);
nz_col = G(:,6);

%Pulse Durations (dimensionalised, in seconds)
tau_col = G(:,2);

%Rabi Rates (dimensionalised)
Omega_col = G(:,1);

%__________________________________________________________________________

%           ITERATIVELY COMPUTE THE CUMULATIVE OPERATOR, Q, 
%                       FOR PULSE SEGMENT INDEX "k"
%__________________________________________________________________________

    for l=1:k;               %Index time intervals by l
            
        %Calculate the evolution operator for the completed pulse 
        %over the lth interval [tlm1, tl]  
        sig_dot_nhat = nx_col(l)*sigx+ny_col(l)*sigy+nz_col(l)*sigz;
    
        RotationExponent = 0.5*Omega_col(l)*tau_col(l)*sig_dot_nhat;

        P=expm(-ii*RotationExponent);
    
        Q=P*Q;                %Q is the cumulative evolution operator. Since
                              %it is essentially a computational tool, and 
                              %we don't care to store it, it can be defined
                              %in a floating manner, where Q in the lth time
                              %interval is obtained by multiplying it's value
                              %in the (l-1)th time interval by the evolution 
                              %operator describing the evolution under the 
                              %control Hamiltonian during the lth interval 
                              %[tlm1, tl].   
    end

%__________________________________________________________________________

%         CONSTRUCT CONTROL HISTORY MATRIX (lambda): index "k"
%__________________________________________________________________________
 
    %matrix lambda (1st row)
    lambda(1,1)=(1/2)*trace(Q'*sigx*Q*sigx);
    lambda(1,2)=(1/2)*trace(Q'*sigx*Q*sigy);
    lambda(1,3)=(1/2)*trace(Q'*sigx*Q*sigz);
    %matrix lambda (2st row)
    lambda(2,1)=(1/2)*trace(Q'*sigy*Q*sigx);
    lambda(2,2)=(1/2)*trace(Q'*sigy*Q*sigy);
    lambda(2,3)=(1/2)*trace(Q'*sigy*Q*sigz);
    %matrix lambda (3rd row)
    lambda(3,1)=(1/2)*trace(Q'*sigz*Q*sigx);
    lambda(3,2)=(1/2)*trace(Q'*sigz*Q*sigy);
    lambda(3,3)=(1/2)*trace(Q'*sigz*Q*sigz);
    
%NOTE: if k ==0, (i.e. a 1-segment control sequence), then the for loop
%running through l = 1:k is bypassed, so Q is simply the identity, as
%defined in the "COMPUTATIONAL OBJECTS" section. Hence, for a 1-segment
%control sequence, the History Matrix, lambda, is trivially the identity. 
   
 HISTORY_MATRIX_GAMMA_MATRIX = lambda;
end



