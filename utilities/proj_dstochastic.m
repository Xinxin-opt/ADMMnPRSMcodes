function   x = proj_dstochastic(x0,tol,n,P,pinvP,Q,pinvQ,en)
%% min ||X-X0||^2_F s.t.Xe =e, X^e = e, X >=0
%%% equivalently, min ||X-X0||^2_F s.t. X = Y, Z = Y, Xe=e, Ye=e, Z>=0
%%% two block ADMM, treating [X;Z] as a variable and Y as another variable
%%% input: x0 --- initial point
%%%        tol --- tolerance
%%%        n  --- dimension
%%%        P ---- P = kron(enT,In); pinvP = P'/PPT;
%%%        Q --- Q = kron(In,enT); pinvQ = Q'/QQT;
%%%        en --- all one vector
%%% output: x -- x=vec(X), X is a doubly stochastic matrix

beta = 1;
maxiter = 10000;

%if max(x0)<1e-3
%    disp('The input matrix is negative!');
%end
y = x0;
lambda1 = zeros(n^2,1);
lambda2 = zeros(n^2,1);

iter = 0; x=abs(y)+1; z=0;
while  max(norm(x-y),norm(y-z))>tol && iter <= maxiter
    iter = iter+1;
    temp = (x0 + beta*(y+lambda1/beta))/(1+beta);
    x = temp - pinvP*(P*temp-en); %%% x subproblem
    z = max(0,y+lambda2/beta); %%%% z subproblem
    lambda1 = lambda1-0.9*beta*(x - y);
    lambda2 = lambda2-0.9*beta*(z - y);
    temp = (x- lambda1/beta + z - lambda2/beta)/2;
    y = temp(:)-pinvQ*(Q*temp-en); %%%% y subproblem
    lambda1 = lambda1-0.9*beta*(x - y);
    lambda2 = lambda2-0.9*beta*(z - y);
    %if max(norm(x-y),norm(y-z))<tol
    %fprintf('projection onto doubly stochastic matrix stops at iter = %i \n',iter);
    %    break;
    %end
end
%X = reshape(x,n,n);
end   % end of function
