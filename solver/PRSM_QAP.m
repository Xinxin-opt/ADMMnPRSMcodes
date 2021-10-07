function [Y, Out] = PRSM_QAP(A, B, Vhat, Jbar, opts, params)
%function [Y, Out] = PRSM_QAP(A, B, Vhat, Jbar, opts, params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Peaceman-Rachford Splitting Method (PRSM) for the DNN relaxation of QAP
% (QAP) min trace AXBX', s.t. X permutation
% (DNN) min trace(LY),     L = blkdiag(1,kron(B,A))
%       s.t.  Vhat R Vhat' = Y,
%             R psd, trace(R) = n + 1
%             0 <= Y <=1;
%             Y(Jbar) = 0, and Y(1,1) = 1;
%             arrow positions of Y holds the doubly stochasticity  (compact notation?)
%%% max-min of augmented Lagrangian
%%% max_Z min_{R,Y} trace LY  + trace (Z*(Y-VRV')) + beta/2 * ||Y-VRV'||^2
%%% Input:
%           A,B : flows, distance matrices, assumed nonnegative integer
%           Vhat: a matrix with orthonormal columns for the facial reduction
%           Jbar: index set in the gangster constraint set for Yij = 0
%           opts.:  structure: options
%           params.: structure: parameters (beta, gamma, initial points)
%%% Output:
%           Y   :  opt solution
%           Out.:  structure(various outputs) best lower bound; best upper
%                  bound; objective function value; primal residual; dual
%                  residual; etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Last version: Sep 20, 2021

%% parameter setting
if isfield(params,'beta'),     beta = params.beta;        else, beta = 50;      end
if isfield(params,'gamma'),    gamma = params.gamma;      else, gamma = 0.9;      end
%% option setting
if isfield(opts,'maxit'),    maxit = opts.maxit;      else, maxit = 40000;    end
if isfield(opts,'tol'),      tol = opts.tol;          else, tol = 1e-4;     end
if isfield(opts,'verbose'),  verbose = opts.verbose;  else, verbose = true; end
if isfield(opts,'dstol'),    dstol = opts.dstol;    else, dstol = 1e-4; end
if isfield(opts,'eval_ubd'),   eval_ubd = opts.eval_ubd;    else, eval_ubd = false; end
if isfield(opts,'eval_lbd'),   eval_lbd = opts.eval_lbd;    else, eval_lbd = false; end
if isfield(opts,'kkt'),      opts.kkt = opts.kkt;     else, opts.kkt = 0; end
if isfield(opts,'lbdstop'),  opts.lbdstop = opts.lbdstop;     else, opts.lbdstop = -inf; end

%% set sizes
n = length(A);
n2 = n^2;
n2p = n2 + 1;

%% coefficients for upper bound function cal_ubd; only done once here
CoefA = [kron(eye(n),ones(1,n)); kron(ones(1,n),eye(n))];
CoefAd = CoefA(1:end-1,:);
rhsb = ones(2*n,1);
rhsbd = rhsb(1:end-1,:);
linopts = optimoptions('linprog','Algorithm','dual-simplex','Display','off');

%% read L for objective trace(L*Y)
L = blkdiag(0,kron(B,A));   % original obj. for DNN relaxation
Lorig = L; % for evaluating obj
if ~issymmetric(L)
    fprintf('BEWARE: The matrix L is not symmetric! so symmetrizing here');
    L = (L+L')/2;
end

%% preprocessing L : project, shift, scale L %%%%%%%%%%%%%%%%%%%%
%%% step 1/3. project L to get equivalent obj. L  since Y = ProjV*Y*ProjV = VRV'
fVhat = full(Vhat);
L = (fVhat*((fVhat'*L)*fVhat))*fVhat';   %L <-- ProjV L ProjV
L = (L + L')/2;  % does not change obj values

%%% step 2/3. shift L    shiftL is sigma_L in paper
shiftL = max(0,-floor(min(eig(L)))) + 10*n;  % ensure suff. pos def
L = L + shiftL*speye(n2p);

%%% step 3/3. scale matrix L for stability  --- work with scaled problem!
nrmL = norm(L,'fro');
scale = nrmL/n2;     % is alpha/n2   in paper
scale = ceil(scale+0.1); %integer look better
L = L/scale;

%%% to recover correct lbd,  use: scale*lbd - (n+1)*shift
%%%%%%%%%%%%%%%%%%%% end of preprocessing L %%%%%%%%%%%%%%%%%%%%


%% initial for R,Y,Z
Y0 = params.Y0;
Z0 = params.Z0;
Y = Y0; Z = Z0;
nrm_pR = 1;
nrm_dR = 1;
VRV = Y;

%% indices for projecting the arrow position elements of Z using L
%%% so that (L+Z)_ij = 0  for all ij in redundant constraints
Zdc = logical([false        true(1,n2)
    true(n2,1)     speye(n2) ]);

%%% vectors to use for updating Z
dcL = -L(Zdc);

Z(Zdc) = dcL;  % first dual iterate

%% Arrays for saving results
feas = zeros(maxit,1);
obj = zeros(maxit,1);
hist_pR = zeros(maxit,1);
hist_dR = zeros(maxit,1);


%% initial results
nstall = 0;   % norm error measures do not change
ustall = 0;   % upper bound measures do not change
lstall = 0;   % lower bound measures do not change
cntubd = 0;   % counts number of times ubd, upper bnd, called cal_ubd
ubest = trace(A*B);   % for X == I ensure upper bound
lbest = 0;    % a lower bound
lbestu = 0;   % for unrounded lower bound
iter = 0;
bestiter = 0; % save the iter where ubest, lbest changed last

Out.lbd = [lbest];
Out.lbdu = [lbest];   % unrounded lower bounds
Out.ubd = [ubest];

nrm_dZ = 0;
relgap = 2*(ubest-lbest)/(ubest+lbest+1);



%% main algorithm
tic_lbdstop = tic;  % for measuring lbd time to compare with SDPNAL
lbdstoptime = - 100; % if this value is never assigned, it remains negative
ubdtime = 0; %save ubd computation time

if opts.kkt==1
    fprintf('\t iter# \t\t nrm_pR  \t  Opt_R  \t   Opt_Y \t   nrm_dR \t   nrm_dZ \t relgap \t obj L''Y  \t time_toc \n');
    fprintf('\t _____\t\t__________\t__________\t__________\t__________\t__________\t_________\t_________\t_________   \n')
else
    fprintf('\t iter# \t\t nrm_pR  \t   nrm_dR \t   nrm_dZ \t relgap \t obj L''Y  \t time_toc \n');
    fprintf('\t _____\t\t__________\t__________\t__________\t_________\t_________\t_________   \n')
end

%%%%%% stopping conditions:
% nstall increases:  if nrm_pR < tol && nrm_dR < tol  % residual nrms
% ustall increases:  if Out.ubd(end) < ubest    % upper bnd not changing
% lstall increases:  if Out.lbd(end) > lbestu   % lower bnd not changing
% ubest == lbest   % upper and lower bnds the same
lstop = 100;
ustop = 100;
nstop = 100;


start_time = tic;   % start time
%%% data used for proj_dstoch %%%%%%%%%%%%%%%%%%%%%%%
en = ones(n,1);
enT = ones(1,n);
In = eye(n);
P = kron(enT,In);
Q = kron(In,enT);
PPT = P*P';
QQT = Q*Q';
pinvP = P'/PPT;
pinvQ = Q'/QQT;
yD = zeros(n);  % for simplex step for arrow
%%% end data used for proj_dstoch %%%%%%%%%%%%%%%%%%%%
kkt = 0;  % a value in [0,1]; [Yoptimality + Roptimality + (Y-VRV)]/3
flaglbdstophit = false; % flag; hit the lbdstop value before or not
start_iter_time = tic;

while      ((ustall < ustop) || (lstall < lstop)) ... % both u/l bnds stall
        &&  (nstall < nstop) ...   % stop if residual nrms not improving
        &&  (ubest ~= lbest) ...   % stop if zero gap
        &&  (iter < maxit)   && (kkt~=1)
    
    iter = iter + 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% step/STEP 1 of 4 : update/project R using eig and simplex
    %%                    The projection is onto psd with trace=n+1.
    fW = (Y + (Z/beta));
    WVhat = fVhat'*(fW*fVhat);
    WVhat = (WVhat'+WVhat)/2;
    [U,s] = eig(WVhat,'vector');
    Rproj = simplex_proj(s,n+1);  % projection onto simplex
    id = find(Rproj>0);
    sn = length(id);
    if ~isempty(id)
        tempid = fVhat*U(:,id);
        VRV = tempid*spdiags(Rproj(id),0,sn,sn)*tempid'; % Vhat*R*Vhat', R NOT updated
    else
        VRV = zeros(length(fVhat));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% step/STEP 2 of 4 : first update for Z
    Z = Z + gamma*beta*(Y - VRV);
    Z(Zdc) = dcL;   % converges to diag and colrow-1 of -L(2:end)
    Z = (Z+Z')/2;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% step/STEP 3 of 4 : update Y
    %%             : project onto gangster, 0,1, AND blk-diag = Is
    Y = VRV - (L + Z)/beta;
    
    %% projection onto polyhedral and gnagster constraints
    Y(Y>1) = 1;
    Y(Y<0) = 0;
    Y(Jbar) = 0; Y(1,1) = 1; %gangster
    
    %% projection onto arrow constraint and doubly stochastic matrices
    yD(:) = (VRV((n2p+2:n2p+1:end)') + 2*VRV(2:end,1))/3; % row vector with avg
    yD(:) = proj_dstochastic(yD(:),dstol,n,P,pinvP,Q,pinvQ,en);
    Y((n2p+2:n2p+1:end)') = yD(:);
    Y(2:end,1) = yD(:);
    Y(1,2:end) = yD(:)';
    
    Y = (Y+Y')/2;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% step/STEP 4 of 4 : update dual variable Z
    pR = Y-VRV;
    Z = Z + gamma*beta*pR;
    Z(Zdc) = dcL;   % converges to diag and colrow-1 of -L(2:end)
    Z = (Z+Z')/2;
    
    %% Supplementary Computations
    if verbose
        nrm_dZ = norm(Z-Z0,'fro');   % measures change after BOTH Z updates
        Z0 = Z;   % for total change in Z after both R,Y updates
    end
    
    dR = Y - Y0; Y0 = Y;
    nrm_pR = norm(pR,'fro');
    feas(iter) = nrm_pR/norm(Y,'fro'); % some relative distance measurement
    
    obj(iter) = Lorig(:)'*Y(:); % objective function value at current iterate
    
    nrm_dR = beta*norm(dR,'fro');
    hist_pR(iter) = nrm_pR;
    hist_dR(iter) = nrm_dR;
    if mod(iter,100) == 0
        relgap = 2*(ubest-lbest)/(ubest+lbest+1);
        if opts.kkt==1  % if KKT condition is a stopping criterion
            %%%%%%%%% optimality check (KKT) -- see the write-up for the formula %%%%%%%%%
            %%%%% optimality condition in Y
            %%%% Y == projection_{set Y} (Y − L − Z)
            LZ = L + Z;
            YLZ  = Y-LZ;    % project this to the set Y
            yD(:) = (YLZ((n2p+2:n2p+1:end)') + 2*YLZ(2:end,1))/3;
            YLZ(YLZ>1) = 1; %proj on polyhedtal set
            YLZ(YLZ<0) = 0;
            YLZ(Jbar) = 0;  %proj on gangster
            YLZ(1,1) =1;
            % project onto DS on the arrow positions
            yD(:) = proj_dstochastic(yD(:),dstol,n,P,pinvP,Q,pinvQ,en);
            YLZ((n2p+2:n2p+1:end)') = yD(:);
            YLZ(2:end,1) = yD(:);
            YLZ(1,2:end) = yD(:)';
            
            optimalityY = norm(YLZ-Y,'fro');   % Y related KKT condition (residual)
            
            %%%%% optimality condition in R
            %%%%% R == projection_{set R}(R + V'ZV)
            R = Vhat'*VRV*Vhat;  R = 0.5*(R+R');
            VZV = Vhat'*Z*Vhat;
            RVZV = R +  VZV;  RVZV = 0.5*(RVZV+RVZV');  % project this to the set R
            [UF,sF] = eig(RVZV,'vector');
            Rproj = simplex_proj(sF,(n+1));
            id = find(Rproj>0);
            sn = length(id);
            tempid = UF(:,id);
            PRVZV = tempid*spdiags(Rproj(id),0,sn,sn)*tempid';
            
            optimalityR = norm(R-PRVZV,'fro') ;  % R related KKT condition (residual)
            
            %%% flag to see if each condition is satisfied
            Rflag = optimalityR<tol; % R condition
            Yflag = optimalityY<tol; % Y condition
            pRflag = nrm_pR<tol;     % pR condition
            
            kkt = (Rflag+Yflag+pRflag)/3; % if kkt=1, then kkt condifion is satisdies
        end  %% end of  'if opts.kkt==1'
        if verbose
            if opts.kkt==1
                fprintf('\t %i \t\t%5.4e \t%5.4e \t%5.4e \t%5.4e \t%5.4e \t%5.4e \t%5.4e \t%5.4e\n',...
                    iter,nrm_pR,optimalityR,optimalityY,nrm_dR,nrm_dZ,relgap,obj(iter),...
                    toc(start_iter_time));
            else   % if KKT condition is not a stopping criterion
                fprintf('\t %i \t\t%5.4e \t%5.4e \t%5.4e \t%5.4e \t%5.4e \t%5.4e\n',...
                    iter,nrm_pR,nrm_dR,nrm_dZ,relgap,obj(iter),...
                    toc(start_iter_time));
            end  % end of   'if opts.kkt==1'
        end % end of 'if verbose'
    end %   end of   'if mod(iter,100) == 0'
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%   update stopping criteria
    if nrm_pR < 1e-1 && nrm_dR < 1e-1
        %% find lbd, ubd in Out.lbd, Out.ubd
        if eval_ubd && mod(iter,100) == 0
            ubdtimetemp = tic;  % measure the ubd computation time
            if ustall < ustop
                cal_ubd;
                if Out.ubd(end) < ubest
                    ubest = Out.ubd(end);
                    ustall = 0;  % reset stall count
                    lstall = 0;
                    nstall = 0;
                    bestiter = iter;
                else
                    ustall = ustall + 1;
                end
            end
            ubdtime = ubdtime + toc(ubdtimetemp);  % accumulate ubd timd
        end  % end of  'eval_ubd && mod(iter,100) == 0'
        if eval_lbd && mod(iter,100) == 0
            if lstall < lstop
                cal_lbd;
                if Out.lbd(end) > lbestu   % unrounded lower bounds
                    lbest = Out.lbd(end);
                    lbestu = Out.lbd(end);
                    lstall = 0;  % reset stall count
                    ustall = 0;
                    nstall = 0;
                    bestiter = iter;
                else
                    lstall = lstall + 1;
                end
            end
            if (opts.lbdstop <= lbest) && (~flaglbdstophit) % if lbest is just as good as lbdstop
                lbdstoptime = toc(tic_lbdstop);  % measure time
                flaglbdstophit = true;           % turn on the switch so that subsequent iteration does not get in this 'if'
            end
        end  % end of 'if eval_lbd && mod(iter,100) == 0'
        %% checking stopping condition
        if nrm_pR < tol && nrm_dR < tol
            nstall = nstall + 1;
        else
            nstall = 0;  % initialize the count
        end
    end  % end of 'if nrm_pR < 1e-1 && nrm_dR < 1e-1 '
    
end % end of main iteration, while loop
end_iter_time = toc(start_iter_time);
Out.toc = toc(start_time);

if iter>0
    cal_ubd
    ubest = min(Out.ubd);
    cal_lbd
    lbest = max(Out.lbd);
    %%%%%%%%% optimality check
    if opts.kkt    % if KKT opt cond was used for stopping criteria
        %%%% optimality condition in Y (residual)
        %%%% Y == projection_{set Y} (Y − L − Z)
        LZ = L+Z;
        YLZ  = Y-LZ;
        %%% projection on set Y
        yD(:) = (YLZ((n2p+2:n2p+1:end)') + 2*YLZ(2:end,1))/3;
        YLZ(YLZ>1) = 1;  % polyhedral set
        YLZ(YLZ<0) = 0;
        YLZ(Jbar) = 0;   % ganster
        YLZ(1,1) =1;
        % DS projection on the arrow positions
        yD(:) = proj_dstochastic(yD(:),dstol,n,P,pinvP,Q,pinvQ,en);
        YLZ((n2p+2:n2p+1:end)') = yD(:);
        YLZ(2:end,1) = yD(:);
        YLZ(1,2:end) = yD(:)';
        optimalityY = norm(YLZ-Y,'fro');
        
        %%%%% optimality condition in R (residual)
        %%%%% R == projection_{set R}(R + V'ZV)
        R = Vhat'*VRV*Vhat;
        R = 0.5*(R+R');
        VZV = Vhat'*Z*Vhat;
        RVZV = R +  VZV;
        RVZV = 0.5*(RVZV+RVZV');
        [UF,sF] = eig(full(RVZV),'vector');
        Rproj = simplex_proj(sF,(n+1));
        id = find(Rproj>0);
        sn = length(id);
        tempid = UF(:,id);
        PRVZV = tempid*spdiags(Rproj(id),0,sn,sn)*tempid';
        optimalityR = norm(R-PRVZV,'fro') ;
        
        % print the final interate information after the while loop
        fprintf('\t %i \t\t%5.4e \t%5.4e \t%5.4e \t%5.4e \t%5.4e \t%5.4e \t%5.4e \t%5.4e\n',...
            iter,nrm_pR,optimalityR,optimalityY,nrm_dR,nrm_dZ,...
            relgap,obj(iter),end_iter_time);
        fprintf('\nnstall, ustall, lstall indicates stalling for norm, upper lower errors resp.\n')
        fprintf('at end of while: kkt = %1.1f, nstall = %i,ustall = %i,lstall = %i,\n                 iter = %i,[lbest ubest] = [%g %g]\n',...
            kkt,nstall,ustall,lstall,iter,lbest,ubest)
    else
        % print the final interate information after the while loop
        fprintf('\t %i \t\t%5.4e \t%5.4e \t%5.4e \t%5.4e \t%5.4e \t%5.4e\n',...
            iter,nrm_pR,nrm_dR,nrm_dZ,relgap,obj(iter),end_iter_time);
        fprintf('\nnstall, ustall, lstall indicates stalling for norm, upper lower errors resp.\n')
        fprintf('at end of while:  nstall = %i,ustall = %i,lstall = %i,\n                  iter = %i,[lbest ubest] = [%g %g]\n',...
            nstall,ustall,lstall,iter,lbest,ubest)
    end  % end of  'if opts.kkt'
end  % end of  'if iter>0'


%% save output
Out.obj = obj(1:iter);   % trace(L*Y)
Out.iter = iter;         % total number iterations
Out.feas = feas(1:iter); % residual nrm_pR/norm(Y,'fro');
Out.pr = hist_pR(1:iter);% primal residual, norm(Y-VRV,'fro')
Out.dr = hist_dR(1:iter);% dual residual, norm(Y-Y0,'fro')
Out.Z = Z;               % dual variable
Out.R = Vhat'*VRV*Vhat;  % primal variable R
Out.Vhat = Vhat;         % Vhat
Out.bestiter = bestiter; % last iteration that yield best bound
Out.ubest = ubest;       % best upper bound
Out.lbest = lbest;       % best lower bound
Out.L = L;               % modified objective function data, trace(L*Y), after scaling and shifting
Out.Lorig = Lorig;       % original objective function data, trace(Lorig*Y)
Out.scale = scale;       % scaling factor of the objective
Out.shiftL = shiftL;     % shifting parameter of the objective
Out.ubdtime = ubdtime;   % time spent on computing upper bounds
Out.lbdstoptime = lbdstoptime;  % first time when the solver lower bound met the user provided lower bound


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% sub function for calculating the lower bound
    function cal_lbd
        % calculate lower bound
        
        Zp = Z;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%% Y related terms in Lagrangain %%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % note : we skip the DS projection due since arrow(L+Z) = 0
        Yp = zeros(n2+1);
        LZp = L + Zp;   % changed here
        
        Yp(LZp < 0) = 1;
        Yp(Jbar) = 0; Yp(1,1) = 1;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%% R related terms in Lagrangain %%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        VZV = Vhat'*Zp*Vhat; %%%% Zp should be here
        VZV = (VZV+VZV')/2;
        Zeigmax = max(eig(VZV));
        
        
        lbd =  sum(sum(LZp.*Yp)) - (n+1)*Zeigmax;  % lbd for scaled probl.
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Out.lbdu = [Out.lbdu; (scale*lbd - shiftL*(n+1))];
        % subtracting 1e-6 in case of roundoff error
        scaled_lbd = ceil((scale*lbd - shiftL*(n+1))-(1e-6) );
        
        if issymmetric(A) && issymmetric(B)
            if mod(scaled_lbd,2)==1   % if the best lbd is odd, then add 1
                scaled_lbd = scaled_lbd + 1 ;
            end
        end
        
        Out.lbd = [Out.lbd; scaled_lbd ];
        
    end    % end of  'cal_lbd'

%% sub function for calculating the upper bound
    function cal_ubd
        % calculate upper bound
        % There are three approaches taken. They all use projection onto
        % the set of permutation matrices from a special point.
        cntubd = cntubd+1;  % number of calls of cal_ubd
        
        Ys = (Y'+Y)/2;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %(approach 1) get feasible solution with first col of Y
        Xloc_hat = full(reshape(Ys(2:end,1),n,n));  % use first col of original Y11
        [xloc,~] = linprog(-Xloc_hat(:),[],[],...
            CoefAd,rhsbd,zeros(n2,1),ones(n2,1),[],linopts);
        Xloc1 = reshape(xloc,n,n);
        
        feas_obj1 = (trace(A*Xloc1*B*Xloc1'));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %(approach 2) get feasible solution with first dominant eigenvector of Y
        
        [u,db] = eig(full(Ys));
        db = diag(db);
        db = db(end:-1:1);
        u  = u(:,end:-1:1);
        Yloc_hat = db(1)*(u(:,1)*u(:,1)');  %u(:,1) is the first dominant eigenvector
        Xloc_hat = reshape(Yloc_hat(2:end,1),n,n);
        
        [xloc,~] = linprog(-Xloc_hat(:),[],[],...
            CoefAd,rhsbd,zeros(n2,1),ones(n2,1),[],linopts);
        Xloc2 = reshape(xloc,n,n);
        
        feas_obj2 = (trace(A*Xloc2*B*Xloc2'));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %(approach 3) perturbation approach
        feas_obj = inf;
        feas_obj3 = inf;
        
        numrank = sum( db > .1);  % determine how much we would consider as numerical rank
        
        if numrank>1
            U = u(:,1:numrank);   % we recycle u obtained above
            rtemp = max(1, min(3*ceil(log(n)),ubest-lbest));
            feas_obj3_save = zeros(rtemp,1);
            
            for jj = 1:rtemp
                rv = rand(numrank,1);
                rv = db(1:numrank).*sort(rv,'descend');
                
                rv = rv/norm(rv);
                
                Yv = U*rv;  % sample vector
                Yv = Yv(2:end);
                [xloc3,~] = linprog(-Yv(:),[],[],...
                    CoefAd,rhsbd,zeros(n2,1),ones(n2,1),[],linopts);
                Xloc3 = reshape(xloc3,n,n);
                feas_obj3  =(trace(A*Xloc3*B*Xloc3'));
                feas_obj3_save(jj) = feas_obj3;
                if feas_obj > feas_obj3
                    feas_obj = feas_obj3;
                end
                
            end
            
            feas_obj3 = feas_obj;
            
        end %%% end if numrank>1
        
        [minfeas,~] = min([feas_obj1,feas_obj2,feas_obj3]);
        
        Out.ubd = [Out.ubd;floor(minfeas+0.01 )];
        
    end %end of  'cal_ubd'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end   % of main function  PRSM_QAP_KKT

