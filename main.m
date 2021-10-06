%%% file is main.m
%%%     runs  qrun_tests.m to prepare to call the solver
%%%% ADMM or PRSM chosen below
%%%% opts/params set here


%% Initialize
clear
addpath(genpath('.'));  % adds path for utilities/solvers/...

onprofile = false;   % profile/debug the run true or false
if onprofile
    profile clear
    profile on
end

%seed = 100;  % since randomization is used in upper bound
%saverng = rng(seed);
saverng = rng('shuffle');


%% run the user's own problem; random problem, or chosen problem from data set
%% datatype = 1,2,3
%% If the user loads their own problem,
%% the data is A,B POS. INTEGER valued, n by n symmetric matrices



%%%%%%%%%%%options follow, choose solver PRSM or ADMM
opts = [];
opts.maxit = 40000;    % upper bound on iterations; defualt = 40000
opts.tol = 1e-4;       % tolerance; default = 1e-4
opts.dstol = 1e-4;     % tolerance for doubly-stochastic set projection; PRSM only; default = 1e-4
opts.verbose = true;   % true=displaying all iteration progress; default = true
opts.PRSM = true;      % run PRSM algorithm/solver or use ADMM with = false; default = true
opts.kkt = false;      % include kkt conditions as stopping criterion; default = false (suggested: opts.kkt=false for small, opts.kkt=true for medium and large)
opts.eval_ubd = true;  % calculate upper and upper bounds; default = true
opts.eval_lbd = true;  % calculate upper and lower bounds; default = true
opts.lbdstop = 0;      % a known lower bound from the user; the function records the first time when opts.lbdstop<=lbest occurs; PRSM option only; ; default = -inf
params = [];


datatype = 1;   % 1,2, or 3
% datatype 1 = user provided instance, 2 = random instance, 3 = QAPLIB
datafilename = 'userdataname'; % current mat file available
filename = 'results';  % filename for saving all objects from qrun_tests.m; data will be saved at results folder

if datatype == 1     %%%  run the user's own problem
    load(datafilename)  % contains symmetric A,B and dim n
    save('datafilename','A','B','n')
    probnumbers = 100;  % used to distinguish with qaplib numbers
    fprintf('\nNew user problem; calling qrun_tests with:')
elseif datatype == 2 %%% generate a random problem
    n = 12;     % pick the size
    A = triu(round(10*rand(n)),1); A=A+A';
    B = triu(round(10*rand(n)),1); B=B+B';
    save('datafilename','A','B','n')
    probnumbers = 100;
    fprintf('\nNew random problem; calling qrun_tests with:')
else                 %%% run a chosen problem from data set
    fprintf('\nNew qaplib problem; calling qrun_tests with:')
    probnumbers = [15:16];  % problem number runs from 1 to 84
end
fprintf('   max # iterations = %i; tolerance = %.2e \n',opts.maxit,opts.tol)



%% Run tests
[Y,Out] = qrun_tests(probnumbers, filename, opts, params);

if onprofile
    profile report
end
