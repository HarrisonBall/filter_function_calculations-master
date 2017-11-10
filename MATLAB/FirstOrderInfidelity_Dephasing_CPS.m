
function [FIRST_ORDER_DEPHASING_INFIDELITY_CPS] = FirstOrderInfidelity_Dephasing_CPS(CPS_matrix,tau,alpha,fc,f0,p)

% ---------------------------------------------------------------------
%                               TEST
% ---------------------------------------------------------------------

% clear all
% alpha =20;
% fc = 200;
% f0 = 4;
% p = 0;
% tau = 1000E-6;
% CPS_matrix = CPS_DCG(tau);
% FirstOrderInfidelity_Dephasing_CPS(CPS_matrix,tau,alpha,fc,f0,p)

%_________________________________________________________________________
        
%                            FUNCTION DESCRIPTION
%_________________________________________________________________________

%This function calculates the first order infidelity of a qubit interacting
%with a dephasing noise field. The first order infidelity is calculated 
%according to the first term in Eq. 4 of [arXiv:1410.1624]. Namely, as a 
%(two-sided) overlap integal between the dephasing PSD, S_z(w) and the 
%dephasing-quadrature filter function F_z(w). 

%The (two-sided) dephasing PSD required for the above calculation is assumed 
%to be derived from engineered phase noise on a carrier resonant with the 
%qubit transition, as described in [arXiv:1403.4632]. In this case, the 
%dephasing field PSD is given, up to a constant, by the first time derivative 
%of the engineered phase noise field. The required constant is picked up 
%by (1) moving to an interaction picture which enables us to recast the 
%Hamiltonian of the system in the form implicitely assumed when computing 
%the infidelity via the overlap integral as in Eq. 4 of [arXiv:1410.1624], 
%and (2) by keeping track of the definition of 1-sided vs 2-sided PSDs. 

%The correct scaling factor is detailed in the notes 
%[Ball2016] Scaling factors for engineered dephasing noise. 

%The inputs to this function are:

%   fc          cutoff frequency of phase-noise field
%   f0          increment between frequency components in phase-noise field
%   alpha       power scaling factor
%   p           power law parameter 
%   tau         control sequence duration
%   CPS_matrix  CPS matrix for control sequence

%The noise PSD takes the form of a discrete frequency comb, with teeth 
%separated by f0, with highest-frequency component at frequency fc, 
%and with power law specified by p:

%                -------------------------------------
%                  NOTE: NOISE POWER LAW  (see notes)
%                -------------------------------------

%                 ________________________________________________ 
    
%     noise      |   1/f^2  |    1/f    |   white   |   Ohmic   |
%     type        ________________________________________________ 

%     p          |    -2    |    -1     |     0     |      1    |
%                 ________________________________________________  

%____________________________________________________________________
        
%                   RESCALE TO ANGULAR FREQUENCY UNITS 
%____________________________________________________________________

%Comb Spacing 
w0 = 2*pi*f0;

%PSD Cutoff Frequency
wc = 2*pi*fc;

%____________________________________________________________________
        
%                 PHASE-TO-DEPHASING SCALING FACTOR 
%____________________________________________________________________


A = 1/4;
B = 1/(2*pi);
C = 2;
D = pi*alpha^2*w0^2/2;

Infidelity_O1_ScalingFactor = A*B*C*D;

%____________________________________________________________________
        
%    DISCRETIZE (UNSCALED) NOISE SPECTRAL COMPONENTS UP TO CUTOFF   
%____________________________________________________________________

%Vector of discrete noise spectral components up to discretized cutoff
J = ceil(wc/w0);
J_vec = [1:J]';
NoisePowerLaw = J_vec.^p;
w_j = w0*J_vec;

%____________________________________________________________________
        
%      CONSTRUCT VECTOR OF FILTER FUNCTION EVALUATIONS @ EACH w_j   
%____________________________________________________________________

FF_comb = FF_CPS_Dephasing(CPS_matrix,w_j,tau);

%____________________________________________________________________

%       CALCULATE THE FIRST ORDER INFIDELITY BY COMPUTING (DISCRETE)
%               OVERLAP INTEGRAL (i.e. SUM) BETWEEN 
%            FILTER FUNCTION AND PRESCRIBED NOISE PSD. 
%____________________________________________________________________

InfidelitySummandFactor = (NoisePowerLaw.*FF_comb)./(w_j.^2); 
Infidelity = Infidelity_O1_ScalingFactor*sum(InfidelitySummandFactor);
Chi = 2*Infidelity;
   
%____________________________________________________________________

%                        OUTPUT GATE INFIDELITY 
%____________________________________________________________________

FIRST_ORDER_DEPHASING_INFIDELITY_CPS = Infidelity;
end





