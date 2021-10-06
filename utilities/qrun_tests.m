function [Y11,Out11] = qrun_tests(probnumbers, filename, opts, params)
%function  qrun_tests(probnumbers, filename, opts, params)
%This function is called from main.m file
%This sets up options, calls/loads the correct problem
%calls the solver PRSM or ADMM; and outputs the final table in latex format.
% input:
%         probnumbers --- problem number in QAP data set
%         filename --- output tex file name
%         opts --- options, opts.maxit;opts.tol; opts.cal_bd; opts.verbose;
%                  opts.PRSM; opts.kkt; etc
%         params --- params.gamma; params.beta; params.Y0; params.Z0;
% output:
%         Y11 --- output variable Y from DNN relaxation
%         Out11--- structure(various outputs)


%% Set seed and choose test data set
if isfield(opts,'verbose'), verbose = opts.verbose;
else, opts.verbose = true; verbose = true; end
if isfield(opts,'PRSM'), PRSM = opts.PRSM; else, PRSM = true; end
if ~isfield(params,'gamma'), params.gamma = .9;  end

numProbs = length(probnumbers);
bestlbds = zeros(numProbs,1);
bestubds = zeros(numProbs,1);
relgaps = zeros(numProbs,1);
optgaps = zeros(numProbs,1);
iters = zeros(numProbs,1);
timessec = zeros(numProbs,1);



%% Choose test data set
Optimal_values;   % load KNOWN optimal values from local file; with file names
rownumb = 0;
%%%%%%%%%%%%%%% data is available in folder data;
for i = probnumbers
    rownumb = rownumb+1;  % row number in table
    % small size
    %%% start chr instance
    %%% One matrix A is the adjacency matrix of a weighted tree the other B that
    %%% of a complete graph.  A === F  and  B === D
    %%% off-diagonal elements/distances in D can be zero
    %%% https://coral.ise.lehigh.edu/data-sets/qaplib/qaplib-problem-instances-and-solutions/#CB
    if     i==1;     loadfile='chr12a';      elseif i==2;    loadfile='chr12b';
    elseif i==3;     loadfile='chr12c';      elseif i==4;    loadfile='chr15a';
    elseif i==5;     loadfile='chr15b';      elseif i==6;    loadfile='chr15c';
    elseif i==7;     loadfile='chr18a';      elseif i==8;    loadfile='chr18b';
    elseif i==9;     loadfile='chr20a';      elseif i==10;    loadfile='chr20b';
    elseif i==11;    loadfile='chr20c';
        %%% end chr instance
        %%% start els instance
        %%A.N. Elshafei [Elshafei:77]
        %%The data describe the distances of 19 different facilities of a
        %%hospital and the flow of patients between those locations. The
        %%optimal solution was first found by [Mautor:92].
        %%  A === D  and  B === F
    elseif i==12;    loadfile='els19';
        %%% end els instance
        %%% start esc instance
    elseif i==13;    loadfile='esc16a';      elseif i==14;    loadfile='esc16b';
    elseif i==15;    loadfile='esc16c';      elseif i==16;    loadfile='esc16d';
    elseif i==17;    loadfile='esc16e';
    elseif i==18;    loadfile='esc16f'; fprintf('WARNING: skip #18!!\n');
    elseif i==19;    loadfile='esc16g';      elseif i==20;    loadfile='esc16h';
    elseif i==21;    loadfile='esc16i';      elseif i==22;    loadfile='esc16j';
        %%% end esc instance
    elseif i==23;    loadfile='had12';       elseif i==24;    loadfile='had14';
    elseif i==25;    loadfile='had16';       elseif i==26;    loadfile='had18';
    elseif i==27;    loadfile='had20';       elseif i==28;    loadfile='nug12';
    elseif i==29;    loadfile='nug14';       elseif i==30;    loadfile='nug15';
    elseif i==31;    loadfile='nug16a';      elseif i==32;    loadfile='nug16b';
    elseif i==33;    loadfile='nug17';       elseif i==34;    loadfile='nug18';
    elseif i==35;    loadfile='nug20';       elseif i==36;    loadfile='rou12';
    elseif i==37;    loadfile='rou15';       elseif i==38;    loadfile='rou20';
    elseif i==39;    loadfile='scr12';       elseif i==40;    loadfile='scr15';
    elseif i==41;    loadfile='scr20';       elseif i==42;    loadfile='tai10a';
    elseif i==43;    loadfile='tai12a';      elseif i==44;    loadfile='tai15a';
    elseif i==45;    loadfile='tai17a';      elseif i==46;    loadfile='tai20a';
        % medium size
    elseif i==47;    loadfile='chr22a';      elseif i==48;    loadfile='chr22b';
    elseif i==49;    loadfile='chr25a';      elseif i==50;    loadfile='esc32a';
    elseif i==51;    loadfile='esc32b';      elseif i==52;    loadfile='esc32c';
    elseif i==53;    loadfile='esc32d';      elseif i==54;    loadfile='esc32e';
    elseif i==55;    loadfile='esc32g';      elseif i==56;    loadfile='esc32h';
    elseif i==57;    loadfile='kra30a';      elseif i==58;    loadfile='kra30b';
    elseif i==59;    loadfile='kra32';      elseif i==60;    loadfile='nug21';
    elseif i==61;    loadfile='nug22';      elseif i==62;    loadfile='nug24';
    elseif i==63;    loadfile='nug25';      elseif i==64;    loadfile='nug27';
    elseif i==65;    loadfile='nug28';      elseif i==66;    loadfile='nug30';
    elseif i==67;    loadfile='ste36a';      elseif i==68;    loadfile='ste36b';
    elseif i==69;    loadfile='ste36c';      elseif i==70;    loadfile='tai25a';
    elseif i==71;    loadfile='tai30a';      elseif i==72;    loadfile='tai35a';
    elseif i==73;    loadfile='tai40a';      elseif i==74;    loadfile='tho30';
    elseif i==75;    loadfile='tho40';
        % large size
    elseif i==76;    loadfile='esc64a';      elseif i==77;    loadfile='sko42';
    elseif i==78;    loadfile='sko49';      elseif i==79;    loadfile='sko56';
    elseif i==80;    loadfile='sko64';      elseif i==81;    loadfile='tai50a';
    elseif i==82;    loadfile='tai60a';      elseif i==83;    loadfile='tai64c';
    elseif i==84;    loadfile='wil50';
    elseif i==100;    loadfile='datafilename';
    end  % end of loading problem in .mat
    
    load([loadfile '.mat'],'A','B','n')   % loads  A,B,n
    if PRSM
        fprintf('\n<strong>[Running PRSM_QAP with problem #i= %i  from file %s]</strong> \n',i,loadfile)
    else
        fprintf('\n<strong>[Running ADMM_QAP with problem #i= %i  from file %s]</strong> \n',i,loadfile)
    end
    if ~isfield(params,'beta'), params.beta = n/3;  end
    %% set the D F
    n2   = n^2;
    In = speye(n);
    En = ones(n);
    
    A(eye(n)==1) = 0;
    B(eye(n)==1) = 0;
    
    %% choose the Vhat for facial reduction
    Vchoice = 3;   % 3 sparsest choice .... 1 is QR, 2 is triangular
    start_iter_time_Vhat = tic;
    if Vchoice==1
        if verbose
            fprintf('using Vchoice=1 standard V with qr for Vhat \n')
        end
        % using qr
        Veye = [eye(n-1); -ones(1,n-1)];
        [V,~] = qr(Veye,0); %If we do this now, we do not normalize the columns of Vhat later
    elseif Vchoice==2
        if verbose
            fprintf('using  Vchoice=2 uppertriangular V for Vhat \n')
        end
        % using upper triangular ones
        V0 = ones(n,n-1);
        V0 = triu(V0);
        V0(2:n+1:end) = -(1:n-1);
        DV0 = sqrt((1:1:(n-1))+(1:1:(n-1)).^2);
        V0 = V0/diag(DV0);
        V = V0;
    elseif Vchoice==3
        if verbose
            fprintf('using  Vchoice=3 sparsest V by blocks for Vhat \n')
        end
        % using orthog. blocks
        iblk=1;                            % number of current block
        ii=1;                              % current beginning column of block iblk
        temp=(spalloc(n,n-1,n*(n+1)/2));     % matrix for final V
        sizeblocks=zeros(floor(n/2),1);    % initialize vector
        normalize=zeros(floor(n/2),1);    % initialize vector
        dscale=ones(n-1,1);              % scaling vector
        sizeblocks(iblk)=floor(n/(2^iblk));     % size/#cols of blk iblk
        normalize(iblk)=norm(ones(2^iblk,1));   %  normalization scalars
        while sizeblocks(iblk) >= 1         % room for at least one column
            tspeye=speye(sizeblocks(iblk));
            ncols=sizeblocks(iblk);
            dscale(ii:ii+ncols-1)=normalize(iblk);
            temp(:,ii:ii+ncols-1)=...
                [ kron(tspeye(:,1:ncols),[ones(2^(iblk-1),1);-ones(2^(iblk-1),1)])
                zeros(n-sizeblocks(iblk)*2^iblk,ncols) ];
            ii=ii+ncols;
            iblk=iblk+1;
            sizeblocks(iblk)=floor(n/(2^iblk));
            normalize(iblk)=norm(ones(2^iblk,1));
        end
        itest= (sum(diag(temp'*temp)==0));
        if itest > 0
            ttemp=sparse(null(full([temp ones(n,1)]')));
            temp(:,end-itest+1:end)=ttemp(:,1:itest);
        end
        dd=sqrt(diag(temp'*temp));
        Ddd=spdiags(dd,0,n-1,n-1);
        V=temp/Ddd;
        V(abs(V)<1e-12)=0;  % zero out any tiny elements
    end    % end if for Vchoice
    %%% Vhat is created for FR   Y = Vhat*R*Vhat'
    KVV = sparse(kron(V,V));  % V is orthog basis for e^perp
    Vhat = [sqrt(1/2), sparse(1,(n-1)^2); sqrt(1/2)/n*ones(n2,1), KVV];
    
    
    
    testorth = Vhat'*Vhat; %check that Vhat'*Vhat = I
    errornormVhat = ...
        norm(testorth-speye(length(testorth)),'fro')/ ...
        sqrt(length(testorth));
    if verbose
        fprintf('rel. error in Vhat''Vhat-I = %g; sparsity Vhat =  %g \n',...
            errornormVhat,nnz(Vhat)/numel(Vhat));
    end
    
    time_Vhat = toc(start_iter_time_Vhat); % time used for obtaining Vhat
    
    %% Forming the gangster index set Jbar  - this does NOT include Y(1,1)
    
    Es21 = triu(ones(n),1);
    YJ = sparse(kron(In,triu(En,1)) + kron(Es21,In));
    YJ = blkdiag(0,YJ);
    YJ = YJ + YJ';
    Jbar = sparse(YJ~=0);     % fix 0 elements in barycenter to stay 0
    
    
    
    %% set initial points and parameters
    
    Yhat = ([0, 1/n*ones(1,n2);
        1/n*ones(n2,1), 1/n2*ones(n2)+1/(n2*(n-1))*kron(n*In-En,n*In-En)]);
    Yhat(1,1) = 1;
    Yhat(Jbar) = 0; Yhat(1,1)=1;
    Yhat = sparse(Yhat);
    
    
    R0 = Vhat'*Yhat*Vhat;
    R0 = (R0 + R0')/2;
    R0(abs(R0)<1e-12)=0;
    Y0 = Yhat;
    Z0 = Y0 - Vhat*R0*Vhat';
    Z0 = (Z0 + Z0')/2;
    Z0(abs(Z0)<1e-12)=0;
    
    params.Y0 = Y0; params.Z0 = Z0;
    
    if verbose
        disp('-------------------------------------------------');
        disp('params:');
        disp(params);
    end
    
    %% Call the solver
    fprintf('\n\nCalling the solver:\n')
    if PRSM
        [Y11,Out11] = PRSM_QAP(A, B, Vhat, Jbar, opts, params);
    else
        [Y11,Out11] = ADMM_QAP(A, B, Vhat, Jbar, opts ,params);
    end
    
    %% Display outputs
    
    %%%%%%%%%% without Vhat time as we can construct a database with Vhat
    hms = sec2hms(Out11.toc);
    
    
    bestlbd = max(Out11.lbd); % bounds from the algorithm
    bestubd = min(Out11.ubd);
    
    % further tightening of the bounds as opt.vals should be even integer
    if mod(bestlbd,2)==1   % if the best lbd is odd, then add 1
        bestlbd = bestlbd + 1 ;
    end
    
    relgap = 2*(bestubd-bestlbd)/(abs(bestubd+bestlbd)+1);
    
    %%% gap between true optimal value
    if probnumbers < 100
        optgap = 2*(bestubd-all_opt(i))/(abs(bestubd+all_opt(i))+1);
    end
    bestiter = Out11.bestiter;
    
    fprintf('\nFrom qrun_tests:\n Max.Lbd.Val: %g\n Min.Ubd(Feas.)Val: %g\n #Iters : %i\n Solver CPU: %s \n Rel. QAP Gap: %g \n Last-iter for best gap: %i \n\n',...
        bestlbd,bestubd, Out11.iter,hms,relgap,bestiter);
    
    bestlbds(rownumb) = bestlbd;
    bestubds(rownumb) = bestubd;
    relgaps(rownumb) = relgap;
    if probnumbers < 100
        optgaps(rownumb) = optgap;
    end
    iters(rownumb) = Out11.iter;
    timessec(rownumb) = Out11.toc + time_Vhat;
    
    
    save(['results\', filename,'.mat'])    % save results
    if i < 100 && loadfile ~= name_all(i)  % i == 100 for special cases
        fprintf('ERROR: loadfile name and optfile name contradiction\n')
        keyboard
    end
    
    
    
end  % of for loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


