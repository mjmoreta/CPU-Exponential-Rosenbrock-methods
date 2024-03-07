% Program that implements exponential Rosenbrock method given by tableau
%
% 0      |   
% 1/2    |  1/2 phi_{1,2}
%  1     |  phi_1
% -------------------------------------------------------------------------
%        |  phi_1-14 phi_3+36 phi_4    16phi_3-48phi4   -2 phi_3+12phi4  
%
% with the standard method of lines 
% u_t=u_xx+f(t,u), x\in [0,1]
% u(t,0)=g0(t), u(t,1)=g1(t) 
% with f(t,u)=u^2+h(t,u) and with exact solution u(t,x)=cos(t+x)
% For the spatial discretization we have used the standard second-order
% difference squeme. 
% This methods appears in 
% V. T. Luan and A. Ostermann, Parallel exponential Rosenbrock methods, 
% Comput. Math. Appl. 71 (2016) 1137-1150.
% Here, we parallelize in order to calculate the intermediate stages

% Some initial values, such as h, matrix A_h0 and some more that are used
% several times along the program
% N is such that h/1/N is the grid diameter in [0,1]

N=1000;
A=-2*diag([ones(N-1,1)])+diag (ones (N-2,1),1)+diag (ones (N-2,1),-1);
h=1/N;
h2=h^2;
A=A./(h^2);
x=[h:h:1-h]';
hdiv=1/h2;

% n is such that the time step size is k=1/n.
n=5;

% We create two parallel pools on cluster.
parpool(2);

% The program runs for 8 different values of k, from k=1/5, in order to calculate the
% error and the order of the method
for ll=1:8
    k=1/n;
    k12=k/2;
    k22=k^2;
    k33=k^3;
    k16=6/k;
    k2212=12/k22;

    tie(1)=k12;
    tie(2)=k;
    Unaux=zeros(N-1,2);

    % U is the approximation to the exact solucion at time t_n. The initial 
    % u(0,x) value is known
    U=zeros(N-1,1);
    for ii=1:N-1   
        U(ii)=cos(x(ii));
    end

    t=0;
    A=sparse(A);
    Ak=k*A;

    Phh=zeros(N-1,1);
    Phh12=zeros(N-1,1);
    Phh1=zeros(N-1,1);
    Vhn=zeros(N-1,1);
    Un2=zeros(N-1,1);
    Un3=zeros(N-1,1);
    vecb=zeros(N-1,5);
    vecb2=zeros(N-1,3);

    % CPU time calculus starts 
    tstart=tic;
    
    % Boundary values that are needed. They are 0 except for ii=1 and ii=N.
    % They are Dirichlet. uu calculates the boundary values of u     
    Chu1=hdiv*uu(t,0);
    ChuN=hdiv*uu(t,1);
    
    % funh is funcion h(t,x)
    for ii=1:N-1
        Phh(ii)=funh(t,x(ii)); 
    end

    % For the local error r=1. For the global one, r=n
    for r=1:n
        % J is the Jacobian
        J=A+2*diag(U);
        J=sparse(J);
        t12=t+k12;
        tk=t+k;
                        
        % Boundary values that are needed. They are 0 except for ii=1 and ii=N.
        % They are Dirichlet. uu calculates the boundary values of u and uut
        % calculates the boundary values of u_t
        Chu121=hdiv*uu(t12,0);
        Chu12N=hdiv*uu(t12,1);  
    
        Chu11=hdiv*uu(tk,0);
        Chu1N=hdiv*uu(tk,1);      
    
        Chut1=hdiv*uut(t,0);
        ChutN=hdiv*uut(t,1);
                   
        % funh is funcion h(t,x) and function funht is h_t(t,x)
        for ii=1:N-1
            Phh12(ii)=funh(t12,x(ii));
            Phh1(ii)=funh(tk,x(ii));
            Vhn(ii)=funht(t,x(ii));        
        end
    
        Vhn(1)=Vhn(1)+Chut1;
        Vhn(N-1)=Vhn(N-1)+ChutN;  
        
        % vecb and vecb2 contain the vectors that are multiplied by the
        % exponential funcions. The first column is zero, the second one is
        % multiplied by phi_1(c_ikA), the third one by phi_2(c_ikA)...

        % Stages
        Au=A*U;
        Uh2masPh=U.^2+Phh;
        vecb2(:,2)=Au+Uh2masPh;
        vecb2(1,2)=vecb2(1,2)+Chu1;
        vecb2(N-1,2)=vecb2(N-1,2)+ChuN;
       
        vecb2(:,3)=Vhn;   
        
        % Stages calculated in parallel
        parfor ll=1:2
            Unaux(:,ll)=U+phipm(tie(ll),J,vecb2,10^(-10),1,1);
        end   
        Un2=Unaux(:,1);
        Un3=Unaux(:,2);

        % Aproximation to the exact solution at time t_{n+1} 
        vecb(:,2)=Au+Uh2masPh;
        vecb(1,2)=vecb(1,2)+Chu1;
        vecb(N-1,2)=vecb(N-1,2)+ChuN;

        vecb(:,3)=Vhn;
    
        Uh2Ph=-U.^2+Phh;
        aux2=Un2.^2+Phh12-2*U.*Un2;
        aux3=Un3.^2+Phh1-2*U.*Un3;
        vecb(:,4)=(-14*Uh2Ph+16*aux2-2*aux3)/k22-k16*Vhn;
        vecb(1,4)=vecb(1,4)+(-14*Chu1+16*Chu121-2*Chu11)/k22;
        vecb(N-1,4)=vecb(N-1,4)+(-14*ChuN+16*Chu12N-2*Chu1N)/k22;  
    
        vecb(:,5)=(36*Uh2Ph-48*aux2+12*aux3)/k33+k2212*Vhn;
        vecb(1,5)=vecb(1,5)+(36*Chu1-48*Chu121+12*Chu11)/k33;
        vecb(N-1,5)=vecb(N-1,5)+(36*ChuN-48*Chu12N+12*Chu1N)/k33;
   
        U=U+phipm(k,J,vecb,10^(-10),1,1); 
   
        t=tk;        

        % We update some values
        Chu1=Chu11;
        ChuN=Chu1N;
        Phh=Phh1;
             
    end

    % CPU time calculus finishes
    telapsed=toc(tstart)

    % Sol contains the exact solution at time T
    sol=zeros(N-1,1);
    for ii=1:N-1
        sol(ii)=cos(x(ii)+t);
    end

    % Error in the infinite norm
    err=norm(sol-U,inf);

    % Order. When ll=1, as there is not a previous error, it can't
    % be calculated. Two consecutive errors are compared
    if ll==1
        err
        err0=err;
    else
        [err log2(err0/err)]
        err0=err;
    end
        
    % The new value of n is 2*n, and the new value of k_n=k/2. It is
    % calculated at the beginning os the next iteration
    n=2*n;

end

delete(gcp('nocreate'));
