=========================================================================
%% solve a doubly nonnegative relaxation for the quadratic assignment problem:  A flows, B distances 
% % %  min tr(AXBX')
% % %  s.t. Xe = e
% % %       X'e = e
% % %       X_{ij} = 0 or 1   (X permutation matrix)
=========================================================================
%% Reference:
% % % A Restricted Dual Peaceman-Rachford Splitting Method for a Strengthed DNN Relaxation for QAP
% % % by Naomi Graham, Hao Hu, Jiyoung Im, Xinxin Li, Henry Wolkowicz
% % % http://www.math.uwaterloo.ca/~hwolkowi/henry/reports/RevisionApr2021newPRSMforQAP.pdf
% % % https://arxiv.org/abs/2006.01529
(as of Sept/21: to appear in INFORMS J. Comput.)

Files and directories in the zip file

README.txt                  this file




main.m.                     To run the code you need to run the script main.m.
This solves one of: the user's own problem; random problem; or chosen problem
from data set. This file calls qrun_tests.m

qrun_tests.m                This sets up options, calls the solver PRSM or ADMM, and outputs relevant information. For explanations on the outputs, see below. 

solver folder:              includes two different codes: 
                            ADMM_QAP.m
                            PRSM_QAP.m

data foler:                 contains instances from QAPLIP,
http://www.mgi.polymtl.ca/anjos/qaplib/inst.html; and
                            Optimal_values.m for best knwon upper bounds from QAPLIP

utilities:                  sec2hms.m, simplex_proj.m, proj_dstochastic.m
                            sec2hms ---  Converts seconds into hours-minimutes-seconds
                            Optimal_values --- contains optimal/best values from each dataset
                            proj_dstochastic --- project a vector onto the set of doubly stochastic matrices
                            simplex_proj --- project a vector onto the simplex

results folder:             folder that stores all objects made in qrun_tests. in .mat file

                            
=========================================================================
qrun_tests.m has two outputs. 
Y:   DNN matrix that is the solver optimal solution. 
Out. : a struct, contains various information
   Out.obj    % history of trace(L*Y) 
   Out.iter   % total number iterations
   Out.feas   % history of residual nrm_pR/norm(Y,'fro') 
   Out.pr     % history of primal residual, norm(Y-VRV,'fro')
   Out.dr     % history of dual residual, norm(Y-Y0,'fro')
   Out.Z = Z;               % final dual variable
   Out.R = Vhat'*VRV*Vhat;  % final primal variable R
   Out.Vhat      % Vhat
   Out.bestiter  % last iteration that yield best bound
   Out.ubest     % best upper bound
   Out.lbest     % best lower bound
   Out.L         % modified objective function data L, trace(L*Y), after scaling and shifting
   Out.Lorig     % original objective function data L, trace(Lorig*Y)
   Out.scale     % scaling factor of the objective
   Out.shift     % shifting parameter of the objective
   Out.ubdtime   % time spent on computing upper bounds
   Out.lbdstoptime  % first time when the solver lower bound met the user provided lower bound (=opts.lbdstop). If is never met, outputs -100.
=========================================================================
