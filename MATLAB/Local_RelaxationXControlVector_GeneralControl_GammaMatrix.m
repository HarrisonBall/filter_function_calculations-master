

function  LOCAL_RELAXATION_X_CONTROL_VECTOR_GENERAL_CONTROL_GAMMA_MATRIX = Local_RelaxationXControlVector_GeneralControl_GammaMatrix(GammaMatrix,l,w_col)

%--------------------------------------------------------------------------
%                                   TEST 
%--------------------------------------------------------------------------
% 
% clear all
% format long
% pitime = 1/4;
% w_col = 0.1;
% GammaMatrix = GammaMatrix_DCG(pitime);
% Local_RelaxationXControlVector_GeneralControl_GammaMatrix(GammaMatrix,l,w_col)

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

%This function outputs the w-domain local control vector associated with
%relaxation noise along the sigma_x axis (i.e. independent of the direction 
%and amplitude of control) for "lth" pulse segment of a piecewise-constant 
%control pulse sequence. That is, the 1st row of the w-domain local control 
%matrix given in Eq. 128 of Ref. [1]). The control Hamiltonian is assumed 
%to satisfy the GENERAL control form. The output of this function is
%copied from the 1st row of the matrix output of the m-file 
%LocalControlMatrix_GeneralControl_GammaMatrix.m. The utility of doing this
%is simply computational speed-up when restricting analysis to noise in the
%sigma_x quadrature. 

%As set out in the supporting mathematica notebooks, GENERAL control
%assumes a control Hamiltonian implementing a sequence of N arbitrary 
%finite-duration rotations of the Bloch vector, corresponding to a sequence 
%of N (arbitary) unitary control operations on a single-qubit. Each pulse
%segment in the sequence is parameterized by a Rabi rate, pulse duration, 
%angle of rotation, and direction of rotation. The direction of rotation 
%is specifed by a unit vector, nhat, pointing in any direction in the Bloch
%sphere:

%           GENERAL CONTROL:        nhat = (nx, ny, nz). 

%That is, each unitary may enact a distinct rotation of the Bloch vector
%about any axis in R^3. A common special case is 

%           EQUATORIAL CONTROL:    nhat = (cos(phi), sin(phi), 0),

%That is, the direction of control is restricted to the equatorial plane 
%of the Bloch sphere (nz = 0), so the direction of the rotation axis may be
%parameterized by a "phase" angle, phi. A further special case is 

%           SINGLE-AXIS (X) CONTROL:   nhat = (1, 0, 0)

%corresponding to the choice phi = 0. 

%The matrix elements for the w-domain local control matrix may be computed 
%by by Fourier transforming the corresponding time-domain local 
%control-matrix elements (see Eq. 128 in Ref. [1])). The expressions for 
%the matrix elements presented in this m-file are taken from the mathematica 
%notebook [H.Ball 2016] local_control_matrix_elements_derivation (MASTER).nb
%(see section "GENERAL CONTROL: w-domain local control-matrix elements stored 
%as functions"). These expressions have been adapted for MATLAB syntax. 

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

%           LOCAL DEPHASING CONTROL VECTOR ELEMENT FORMULAE
%__________________________________________________________________________

G = GammaMatrix;
FRAC = 1./(w_col.^2-G(l,1).^2);

RxxPlwGENERAL= FRAC.*(w_col.^2-G(l,1).^2+G(l,5).^2.*G(l,1).^2+G(l,6).^2.*G(l,1).^2+exp(1i.*w_col.*G(l,2)).*(-w_col.^2.*cos(G(l,3)).*(G(l,5).^2+G(l,6).^2)+1i.*w_col.*sin(G(l,3)).*(G(l,5).^2+G(l,6).^2).*G(l,1)+(-1+G(l,5).^2+G(l,6).^2).*(w_col.^2-G(l,1).^2)));

RxyPlwGENERAL= FRAC.*(-G(l,1).*(1i.*w_col.*G(l,6)+G(l,4).*G(l,5).*G(l,1))+exp(1i.*w_col.*G(l,2)).*(w_col.*G(l,6).*(w_col.*sin(G(l,3))+1i.*cos(G(l,3)).*G(l,1))+G(l,4).*G(l,5).*(w_col.^2.*(-1+cos(G(l,3)))-1i.*w_col.*sin(G(l,3)).*G(l,1)+G(l,1).^2)));

RxzPlwGENERAL= FRAC.*(G(l,1).*(1i.*w_col.*G(l,5)-G(l,4).*G(l,6).*G(l,1))+exp(1i.*w_col.*G(l,2)).*(-w_col.*G(l,5).*(w_col.*sin(G(l,3))+1i.*cos(G(l,3)).*G(l,1))+G(l,4).*G(l,6).*(w_col.^2.*(-1+cos(G(l,3)))-1i.*w_col.*sin(G(l,3)).*G(l,1)+G(l,1).^2)));

%__________________________________________________________________________

%                                OUTPUT
%__________________________________________________________________________

 
LOCAL_RELAXATION_X_CONTROL_VECTOR_GENERAL_CONTROL_GAMMA_MATRIX = [RxxPlwGENERAL RxyPlwGENERAL RxzPlwGENERAL];

end



