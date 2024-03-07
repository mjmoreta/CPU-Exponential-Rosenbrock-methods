% Program that implements method EPIRK4s3 for problem
% u_t=u_xx+f(t,u), x\in [0,1]
% u(t,0)=g0(t), u(t,1)=g1(t) 
% with f(t,u)=u^2+h(t,u) and with exact solution u(t,x)=cos(t+x)
% For the spatial discretization we have used the standard second-order
% difference squeme. 
% This methods has been specially designed for autonomuos problems and used
% special Krylov subroutines KIOPS. The method and the subroutines appears 
% S. Gaudreault, G. Rainwater and M. Tokman, KIOPS: A fast adaptive Krylov 
% subspace solver for exponential integrators. J. Comp. Phys. 372 (2018) 
% 236-255. 

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

% The program runs for 8 different values of k, from k=1/5, in order to calculate the
% error and the order of the method
for ll=1:8
    k=1/n;
    k2=k/2;
    k23=2*k/3;
    k7=7/k;
    k22=k^2;
    k33=k^3;
    k222=2*k22;
    k2218=18/k22;
    tie=[k2,k23,k];
    
    % Variables Y and Yaux contain t and U
    Ynaux=zeros(N,3);
    Y=zeros(N,1);
    
    % U is the approximation to the exact solucion at time t_n. The initial 
    % u(0,x) value is known
    U=zeros(N-1,1);
    for ii=1:N-1    
        U(ii)=cos(x(ii));
    end

    t=0;   
    Y=[t;U];
    t=0;
    A=sparse(A);
    Ak=k*A;

    Phh=zeros(N-1,1);
    Phh2=zeros(N-1,1);
    Phh23=zeros(N-1,1);
    Vhn=zeros(N-1,1);
    Un2=zeros(N-1,1);
    Un3=zeros(N-1,1);
    Yn2=zeros(N,1);
    Yn3=zeros(N,1);
    vecb=zeros(N-1,5);
    vecb2=zeros(N,3);

    % CPU time calculus starts 
    tstart=tic;

    % For the local error r=1. For the global one, r=n
    for r=1:n
        % J is the Jacobian of the non autonomous problem
        J=A+2*diag(U);    
        J=sparse(J);
        t2=t+k2;
        t23=t+k23;  

        % Boundary values that are needed. They are 0 except for ii=1 and ii=N.
        % They are Dirichlet. uu calculates the boundary values of u and uut
        % calculates the boundary values of u_t
        Chu1=hdiv*uu(t,0);
        ChuN=hdiv*uu(t,1);
    
        Chu21=hdiv*uu(t2,0);
        Chu2N=hdiv*uu(t2,1);  
    
        Chu231=hdiv*uu(t23,0);
        Chu23N=hdiv*uu(t23,1);      
    
        Chut1=hdiv*uut(t,0);
        ChutN=hdiv*uut(t,1);
                   
        % funh is funcion h(t,x) and function funht is h_t(t,x)
        for ii=1:N-1
            Phh(ii)=funh(t,x(ii)); 
            Phh2(ii)=funh(t2,x(ii));
            Phh23(ii)=funh(t23,x(ii));
            Vhn(ii)=funht(t,x(ii));        
        end
    
        Vhn(1)=Vhn(1)+Chut1;
        Vhn(N-1)=Vhn(N-1)+ChutN;   

        % Jacobian of the autonomous problem
        Jtilde=[zeros(1,N) ; Vhn J];
        Jtilde=sparse(Jtilde);
        Yh2=Y.^2;
        Uh2=Yh2(2:N);

        % vecb and vecb2 contain the vectors that are multiplied by the
        % exponential funcions. The first column is zero, the second one is
        % multiplied by phi_1(c_ikA), the third one by phi_2(c_ikA)...
        
        % Stages calculus. It is also calculated the part that is
        % multiplied by phi_1 in the numerical approximation to the solution        
        Au=A*U;
        vecb2(1,2)=1;
        vecb2(2:N,2)=Au+Uh2+Phh;
        vecb2(2,2)=vecb2(2,2)+Chu1;
        vecb2(N,2)=vecb2(N,2)+ChuN;
        
        Yaux=kiops(tie,Jtilde,vecb2,10^(-10),1,10,128,true);
        Un2=U+k2*Yaux(2:N,1);
        Un3=U+k23*Yaux(2:N,2);   
 
        aux2=Un2.^2-2*U.*Un2+Phh2;
        aux3=Un3.^2-2*U.*Un3+Phh23;
        
        % Aproximation to the exact solution at time t_{n+1}       
        aux1=-Uh2+Phh;

        vecb(:,4)=(-37*aux1+64*aux2-27*aux3)/k222-k7*Vhn;
        vecb(1,4)=vecb(1,4)+(-37*Chu1+64*Chu21-27*Chu231)/k222;
        vecb(N-1,4)=vecb(N-1,4)+(-37*ChuN+64*Chu2N-27*Chu23N)/k222;  
    
        vecb(:,5)=(63*aux1-144*aux2+81*aux3)/k33+k2218*Vhn;
        vecb(1,5)=vecb(1,5)+(63*Chu1-144*Chu21+81*Chu231)/k33;
        vecb(N-1,5)=vecb(N-1,5)+(63*ChuN-144*Chu2N+81*Chu23N)/k33;
    
        U=U+k*Yaux(2:N,3)+kiops(k,J,vecb,10^(-10),1,10,128,false);

        t=t+k;        
        Y=[t;U];            
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

