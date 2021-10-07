function [Y, Out] = ADMM_QAP(A, B, Vhat, Jbar, opts, params)
%function [Y, Out] = ADMM_QAP(A, B, Vhat, Jbar, opts, params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alternating Direction Method of Mutlipliers (ADMM) for the DNN relaxation
%  of QAP
% (QAP) min trace AXBX',  s.t. X permutation
% (DNN) min trace(LY),   L = blkdiag(1,kron(B,A))
%       s.t.  Vhat R Vhat' = Y,
%             R psd, trace(R) = n + 1
%             0 <= Y <=1;
%             Y(Jbar) = 0, and Y(1,1) = 1;
%%% max-min of augmented Lagrangian
%%% max_Z min_{R,Y} trace LsY  + trace (Z*(Y-VRV')) + beta/2 %norm(Y-VRV')^2
%%% Input:
%%% max_Z min_{R,Y} trace LY  + trace (Z*(Y-VRV')) + beta/2 * ||Y-VRV'||^2
%%% Input:
%           A,B : flows, distance matrices, assumed nonnegative integer
%           Vhat: a matrix with orthonormal columns for the facial reduction
%           Jbar: index set in the gangster constraint set for Yij = 0
%           opts.:  structure: options
%%% Output:
%           Y   :  opt solution
%           Out.:  structure(various outputs) best lower bound; best upper
%                  bound; objective function value; primal residual; dual
%                  residual; etc.
%%% Reference:
% % % D. Oliveira, H. Wolkowicz and Y. Xu, ADMM for the SDP relaxation of the QAP,
% % % Math. Program. Comput.,10(2018),pp. 631-658.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Last version: Sep 20, 2021

%% option setting
if isfield(opts,'maxit'),    maxit = opts.maxit;      else, maxit = 500;    end
if isfield(opts,'tol'),      tol = opts.tol;          else, tol = 1e-1;     end
if isfield(opts,'beta'),     beta = opts.beta;        else, beta = 50;      end
if isfield(opts,'gamma'),    gamma = opts.gamma;      else, gamma = 1;      end
if isfield(opts,'kkt'),      opts.kkt = opts.kkt;     else, opts.kkt = 0; end
if isfield(opts,'eval_ubd'),   eval_ubd = opts.eval_ubd;    else, eval_ubd = false; end
if isfield(opts,'eval_lbd'),   eval_lbd = opts.eval_lbd;    else, eval_lbd = false; end

%% set sizes
n = length(A);
en = ones(n,1);

In = eye(n);

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

%% read L
L = blkdiag(0,kron(B,A));   % original obj. for DNN relaxation
Lorig = L; % for evaluating obj
if ~issymmetric(L)
    fprintf('BEWARE: The matrix L is not symmetric! so symmetrizing here');
    L = (L+L')/2;
    % there are some wrong instances: i=49
end

L = (L + L')/2;  % does not change obj values

%% shift L    shiftL is sigma_L in paper
shiftL = 0; % change   % ADMM did not shift
L = L+shiftL*speye(n2p);

%%%%%%%%%% now scale new L
nrmL = norm(L,'fro');
scale = nrmL/n2;     % is alpha/n2   in paper
scale = ceil(scale+0.1); %integer look better

% scale matrix L for stability  --- work with scaled problem!
L = L/scale;
%%% to recover lbd use:    scale*lbd - (n+1)*shift
%%%%%%%%%% end of scale new L

fVhat = full(Vhat);

%% initial for R,Y,Z
Y0 = params.Y0;
Z0 = params.Z0;
Y = Y0; Z = Z0;   % R is NOT updated
VRV = Y;
nrm_pR=1;
nrm_dR=1;

%% Arrays for saving results
feas = zeros(maxit,1);
obj = zeros(maxit,1);
hist_pR = zeros(maxit,1);
hist_dR = zeros(maxit,1);


%% initial results
nstall = 0;   % norm error measures do not change
ustall = 0;   % upper bound measures do not change
lstall = 0;   % lower bound measures do not change
%ubest = n^2*1000*trace(A*B)/scale;   % for X == I ensure upper bound
ubest = trace(A*B);   % for X == I ensure upper bound
lbest = 0;
lbestu = 0;   % for unrounded lower bound
iter = 0;
bestiter = 0;   % save the iter where ubest, lbest changed last

nrm_dZ = 0;
relgap = 2*(ubest-lbest)/(ubest+lbest+1);

Out.lbd = [lbest];
Out.ubd = [ubest];

%%%%%% stopping conditions:
% nstall increases:  if nrm_pR < tol && nrm_dR < tol  % residual nrms
% ustall increases:  if Out.ubd(end) < ubest    % upper bnd not changing
% lstall increases:  if Out.lbd(end) > lbestu   % lower bnd not changing
% ubest == lbest   % upper and lower bnds the same
lstop = 100;
ustop = 100;
nstop = 100;
kkt = 0;   %a value in [0,1]; [Yoptimality + Roptimality + (Y-VRV)]/3

start_iter_time = tic;

%% main algorithm
start_time = tic;   % start time
ubdtime = 0;  % time used for ubd computations
if opts.kkt==1
    fprintf('\t iter# \t\t nrm_pR  \t  Opt_R  \t   Opt_Y \t   nrm_dR \t   nrm_dZ \t relgap \t obj L''Y  \t time_toc \n');
    fprintf('\t _____\t\t__________\t__________\t__________\t__________\t__________\t_________\t_________\t_________   \n')
else
    fprintf('\t iter# \t\t nrm_pR  \t   nrm_dR \t   nrm_dZ \t relgap \t obj L''Y  \t time_toc \n');
    fprintf('\t _____\t\t__________\t__________\t__________\t_________\t_________\t_________   \n')
end

%%%%%%%%%%main while loop


while (nstall < nstop) && ((ustall < ustop) || (lstall < lstop)) && ...
        (iter < maxit) && (ubest ~= lbest) && (kkt~=1)
    
    iter = iter + 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% step 1 of 3 update R
    fW = (Y + (Z/beta));
    WVhat = fVhat'*(fW*fVhat);
    WVhat = (WVhat'+ WVhat)/2;
    [U,s] = eig(WVhat,'vector');
    %     Rproj = simplex_proj(s,n+1);
    Rproj = s;
    id = find(Rproj>0);
    sn = length(id);
    if ~isempty(id)
        tempid = fVhat*U(:,id);
        VRV = tempid*spdiags(Rproj(id),0,sn,sn)*tempid'; % % Vhat*R*Vhat', R NOT updated
    else
        VRV = zeros(length(fVhat));
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%    step 2 of 3 update Y
    Y = VRV - (L + Z)/beta;
    
    Y(Y>1) = 1;
    Y(Y<0) = 0;
    Y(Jbar) = 0; Y(1,1) = 1; %gangster
    
    Y = (Y+Y')/2;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%    step 3 of 3 update Z
    pR = Y-VRV;
    Z = Z + gamma*beta*pR;
    Z = (Z+Z')/2;
    
    
    %%
    nrm_dZ = norm(Z-Z0,'fro');
    Z0 = Z;
    
    dR = Y - Y0; Y0 = Y;
    nrm_pR = norm(pR,'fro');
    feas(iter) = nrm_pR/norm(Y,'fro'); %some relative distance measurement?
    
    obj(iter) = Lorig(:)'*Y(:); %objective function value at current iterate
    
    nrm_dR = beta*norm(dR,'fro'); %why multiplied by beta? test?????
    hist_pR(iter) = nrm_pR;
    hist_dR(iter) = nrm_dR;
    if mod(iter,100) == 0
        relgap = 2*(ubest-lbest)/(ubest+lbest+1);
        if opts.kkt==1
            
            %%%%%%%%% optimality check (KKT) -- see the write-up for the formula %%%%%%%%%
            %%%%% optimality condition in Y
            %%%% Y == projection_{set Y} (Y − L − Z)
            LZ = L+Z;
            YLZ  = Y-LZ;
            YLZ(YLZ>1) = 1; %proj on polyhedtal set
            YLZ(YLZ<0) = 0;
            YLZ(Jbar) = 0;  %proj on gangster
            YLZ(1,1) =1;
            optimalityY = norm(YLZ-Y,'fro');
            
            %%%%% optimality condition in R
            %%%%% R == projection_{set R}(R + V'ZV)
            R = Vhat'*VRV*Vhat;  R = 0.5*(R+R');
            VZV = Vhat'*Z*Vhat;
            RVZV = R +  VZV;  RVZV = 0.5*(RVZV+RVZV');
            [UF,sF] = eig(RVZV,'vector');
            Rproj = sF;
            id = find(Rproj>0);
            sn = length(id);
            tempid = UF(:,id);
            PRVZV = tempid*spdiags(Rproj(id),0,sn,sn)*tempid';
            optimalityR = norm(R-PRVZV,'fro') ;
            
            %%% flag to see if each condition is satisfied
            Rflag = optimalityR<tol;
            Yflag = optimalityY<tol;
            pRflag = nrm_pR<tol;
            
            kkt = (Rflag+Yflag+pRflag)/3;
            fprintf('\t %i \t\t%5.4e \t%5.4e \t%5.4e \t%5.4e \t%5.4e \t%5.4e \t%5.4e \t%5.4e\n',...
                iter,nrm_pR,optimalityR,optimalityY,nrm_dR,nrm_dZ,relgap,obj(iter),...
                toc(start_iter_time));
        else
            fprintf('\t %i \t\t%5.4e \t%5.4e \t%5.4e \t%5.4e \t%5.4e \t%5.4e\n',...
                iter,nrm_pR,nrm_dR,nrm_dZ,relgap,obj(iter),...
                toc(start_iter_time));
        end  %% end of  'if kkt'
    end %   end of  'mod(iter,100) == 0'
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%   update stopping criteria
    if nrm_pR < 1e-1 && nrm_dR < 1e-1
        if eval_ubd && mod(iter,100) == 0
            %% find lbd, ubd in Out.lbd, Out.ubd
            ubdtimetemp = tic;
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
            ubdtime = ubdtime + toc(ubdtimetemp);
        end  % end of eval_ubd && mod(iter,100) == 0
        
        if eval_lbd && mod(iter,100) == 0
            if lstall < lstop
                cal_lbd;
                if Out.lbd(end) > lbestu   % unrounded bounds
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
        end  % end of eval_lbd && mod(iter,100) == 0
        %% checking stopping condition
        if nrm_pR < tol && nrm_dR < tol
            nstall = nstall + 1;
        else
            nstall = 0;  % initialize the count
        end
    end  %  end of  'nrm_pR < 1e-1 && nrm_dR < 1e-1'
end % end of main iteration, while loop


if iter ==0
    iter = 1;
    Out.lbd = lbest;
    Out.ubd = ubest;
else
    
    cal_ubd
    ubest = min(Out.ubd);
    cal_lbd
    lbest = max(Out.lbd);
    if opts.kkt==1
        %%%% optimality in Y
        LZ = L+Z;
        YLZ  = Y-LZ;
        %%% projection on set Y
        YLZ(YLZ>1) = 1;
        YLZ(YLZ<0) = 0;
        YLZ(Jbar) = 0;
        YLZ(1,1) =1;
        
        optimalityY = norm(YLZ-Y,'fro');
        
        %%%%% optimality in R
        R = Vhat'*VRV*Vhat;
        R = 0.5*(R+R');
        VZV = Vhat'*Z*Vhat;
        RVZV = R +  VZV;
        RVZV = 0.5*(RVZV+RVZV');
        [UF,sF] = eig(full(RVZV),'vector');
        Rproj = sF;
        id = find(Rproj>0);
        sn = length(id);
        tempid = UF(:,id);
        PRVZV = tempid*spdiags(Rproj(id),0,sn,sn)*tempid';
        
        optimalityR = norm(R-PRVZV,'fro') ;
        fprintf('\t %i \t\t%5.4e \t%5.4e \t%5.4e \t%5.4e \t%5.4e \t%5.4e \t%5.4e \t%5.4e\n',...
            iter,nrm_pR,optimalityR,optimalityY,nrm_dR,nrm_dZ,relgap,obj(iter),...
            toc(start_iter_time));fprintf('\nat end of while: kkt = %1.1f, nstall = %i,ustall = %i,lstall = %i,iter = %i,[lbest ubest] = [%g %g]\n',...
            kkt,nstall,ustall,lstall,iter,lbest,ubest)
    else
        fprintf('\t %i \t\t%5.4e \t%5.4e \t%5.4e \t%5.4e \t%5.4e \t%5.4e\n',...
            iter,nrm_pR,nrm_dR,nrm_dZ,relgap,obj(iter),...
            toc(start_iter_time));fprintf('\nat end of while:  nstall = %i,ustall = %i,lstall = %i,iter = %i,[lbest ubest] = [%g %g]\n',...
            nstall,ustall,lstall,iter,lbest,ubest)
    end
end

%% save output
Out.toc = toc(start_time);
Out.obj = obj(1:iter);    % trace(L*Y)
Out.iter = iter;          % total number iterations
Out.feas = feas(1:iter);  % residual nrm_pR/norm(Y,'fro');
Out.pr = hist_pR(1:iter); % primal residual, norm(Y-VRV,'fro')
Out.dr = hist_dR(1:iter); % dual residual, norm(Y-Y0,'fro')
Out.Z = Z;                % dual variable
Out.bestiter = bestiter;  % last iteration that yield best bound
Out.ubest = ubest;        % best upper bound
Out.lbest = lbest;        % best lower bound
Out.L = L;                % modified objective function data, trace(L*Y), after scaling
Out.Lorig = Lorig;        % original objective function data, trace(Lorig*Y)
Out.scale = scale;        % scaling factor of the objective
Out.shiftL = shiftL;      % shifting parameter of the objective
Out.ubdtime = ubdtime;    % time used for ubd computations
Out.lbdstoptime = -100;      % dummy output to work with PRSM output



%% sub function for calculating the lower bound
    function cal_lbd
        % calculate lower bound
        
        That = [-ones(2*n,1),  [kron(In,en'); kron(en',In)] ];
        [Q,~] = qr(That',0);
        Q = Q(:,1:end-1);
        
        Uloc = [Vhat, Q];
        Zloc = Uloc'*Z*Uloc;
        
        W12 = Zloc(1:(n-1)^2+1,(n-1)^2+2:end);
        W22 = Zloc((n-1)^2+2:end,(n-1)^2+2:end);
        W11 = Zloc(1:(n-1)^2+1,1:(n-1)^2+1);
        W11 = (W11+W11')/2;
        [Uw,Dw] = eig(full(W11));
        dw = diag(Dw);
        id = dw<0;
        W11 = Uw(:,id)*diag(dw(id))*Uw(:,id)';
        Zp = Uloc*[W11,W12;W12',W22]*Uloc';
        Zp = (Zp+Zp')/2;
        Yp = zeros(n2+1);
        Yp(L+Zp < 0) = 1;
        Yp(Jbar) = 0; Yp(1,1) = 1;
        Out.lbd = [Out.lbd; ceil(sum(sum((L+Zp).*Yp))*scale)]; %%% unscale
    end   % end of  'cal_lbd'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sub function for calculating the upper bound
    function cal_ubd
        % calculate upper bound
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %(appoach 1) get feasible solution with first dominant eigenvector of Y
        Yloc = (Y'+Y)/2;
        [u,d] = eigs(Yloc,1,'LA');  % we should consolidate the cases to save time on eval computation time
        Yloc_hat = d*(u*u');
        Xloc_hat = reshape(Yloc_hat(2:end,1),n,n);
        
        [xloc,~] = linprog(-Xloc_hat(:),[],[],...
            CoefAd,rhsbd,zeros(n2,1),ones(n2,1),[],linopts);
        Xloc = reshape(xloc,n,n);
        
        feas_obj1 = (trace(A*Xloc*B*Xloc'))/scale;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %(approach 2) get feasible solution with first col of Y
        Xloc = full(reshape(Y(2:end,1),n,n));  % use first col of original Y11
        [xloc,~] = linprog(-Xloc(:),[],[],...
            CoefAd,rhsbd,zeros(n2,1),ones(n2,1),[],linopts);
        Xloc = reshape(xloc,n,n);
        
        feas_obj2 = (trace(A*Xloc*B*Xloc'))/scale;
        
        [minfeas,~] = min([feas_obj1,feas_obj2]);
        
        Out.ubd = [Out.ubd;floor(scale*(minfeas) +0.1 )];
        
        
    end %end of  'cal_ubd'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end   % of main function ADMM_QAP

