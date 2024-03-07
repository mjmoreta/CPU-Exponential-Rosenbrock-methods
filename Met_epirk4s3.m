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
    k8=k/8;
    k9=k/9;
    k34=34/k;
    k22=k^2;
    k33=k^3;
    tie=[k9,k8,k];
    
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
    
    % La variable Y contiene a t y a U
    Y=[t;U];
    
    t=0;
    A=sparse(A);
    Ak=k*A;

    Phh=zeros(N-1,1);
    Phh8=zeros(N-1,1);
    Phh9=zeros(N-1,1);
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
               
        t8=t+k8;
        t9=t+k9;  
            
        % Boundary values that are needed. They are 0 except for ii=1 and ii=N.
        % They are Dirichlet. uu calculates the boundary values of u and uut
        % calculates the boundary values of u_t
        Chu1=hdiv*uu(t,0);
        ChuN=hdiv*uu(t,1);
    
        Chu81=hdiv*uu(t8,0);
        Chu8N=hdiv*uu(t8,1);  
    
        Chu91=hdiv*uu(t9,0);
        Chu9N=hdiv*uu(t9,1);      
    
        Chut1=hdiv*uut(t,0);
        ChutN=hdiv*uut(t,1);
             
        % funh is funcion h(t,x) and function funht is h_t(t,x)
        for ii=1:N-1
            Phh(ii)=funh(t,x(ii)); 
            Phh8(ii)=funh(t8,x(ii));
            Phh9(ii)=funh(t9,x(ii));
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
        Un2=U+k8*Yaux(2:N,2);
        Un3=U+k9*Yaux(2:N,1);   
 
        aux2=Un2.^2-2*U.*Un2+Phh8;
        aux3=Un3.^2-2*U.*Un3+Phh9;
       
        % Aproximation to the exact solution at time t_{n+1}  
        aux1=-Uh2+Phh;
        
        vecb(:,4)=(-434*aux1-1024*aux2+1458*aux3)/k22-k34*Vhn;
        vecb(1,4)=vecb(1,4)+(-434*Chu1-1024*Chu81+1458*Chu91)/k22;
        vecb(N-1,4)=vecb(N-1,4)+(-434*ChuN-1024*Chu8N+1458*Chu9N)/k22;  
    
        vecb(:,5)=(7344*aux1+27648*aux2-34992*aux3)/k33+432*Vhn/k22;
        vecb(1,5)=vecb(1,5)+(7344*Chu1+27648*Chu81-34992*Chu91)/k33;
        vecb(N-1,5)=vecb(N-1,5)+(7344*ChuN+27648*Chu8N-34992*Chu9N)/k33;
    
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


