% Program that implements exponential Rosenbrock method given by tableau
%
% 0      |   
% 3/4    |  3/4 phi_{1,2}
% -------------------------------------------------------------------------
%        |  phi_1-32/9 phi3    32/9 phi_3  
%
% with the standard method of lines 
% u_t=u_xx+f(t,u), x\in [0,1]
% u(t,0)=g0(t), u(t,1)=g1(t) 
% with f(t,u)=u^2+h(t,u) and with exact solution u(t,x)=cos(t+x)
% For the spatial discretization we have used the standard second-order
% difference squeme. 
% This methods appears in 
% V. T. Luan, Fourth-order two-stage explicit exponential integrators for 
% time-dependent PDEs, Appl. Num. Math. 112 (2017), 91-103.

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

% The program runs for 8 different values of k, from k=1/5, in order to 
% calculate the error and the order of the method
for ll=1:8
    k=1/n;
    k34=k*3/4;
    k83=8/(3*k);
    k22=k^2;
    k2329=32/(9*k22);
    k3=k*3;
    
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
    Phh34=zeros(N-1,1);
    Vhn=zeros(N-1,1);
    Un2=zeros(N-1,1);
    Gn1=zeros(N-1,1);
    Gn2=zeros(N-1,1);
    vecb=zeros(N-1,4);
    vecb2=zeros(N-1,3);

    % CPU time calculus starts 
    tstart=tic;
    
    % For the local error r=1. For the global one, r=n
    for r=1:n
        % J is the Jacobian 
        J=A+2*diag(U);
        J=sparse(J);
    
        t34=t+k34;
            
        % Boundary values that are needed. They are 0 except for ii=1 and ii=N.
        % They are Dirichlet. uu calculates the boundary values of u and
        % uut calculates the boundary values of u_t
        Chu1=hdiv*uu(t,0);
        ChuN=hdiv*uu(t,1);
    
        Chu341=hdiv*uu(t34,0);
        Chu34N=hdiv*uu(t34,1);  
    
        Chut1=hdiv*uut(t,0);
        ChutN=hdiv*uut(t,1);
                   
        % funh is funcion h(t,x) and function funht is h_t(t,x)
        for ii=1:N-1
            Phh(ii)=funh(t,x(ii)); 
            Phh34(ii)=funh(t34,x(ii)); 
            Vhn(ii)=funht(t,x(ii));        
        end
    
        Vhn(1)=Vhn(1)+Chut1;
        Vhn(N-1)=Vhn(N-1)+ChutN;   

        Uh2=U.^2;
        
        % vecb2 and vecb contain the vectors that are multiplied by the
        % exponential funcions. The first column is zero, the second one is
        % multiplied by phi_1(c_ikA), the third one by phi_2(c_ikA)...
        
        % Stage     
        vecb2(:,2)=A*U+Uh2+Phh;
        vecb2(1,2)=vecb2(1,2)+Chu1;
        vecb2(N-1,2)=vecb2(N-1,2)+ChuN;
       
        vecb2(:,3)=Vhn;
    
        Un2=U+phipm(k34,J,vecb2,10^(-10),1,1); 
           

        % Aproximation to the exact solution at time t_{n+1}
        vecb(:,2)=vecb2(:,2);   
    
        vecb(:,3)=Vhn;
    
        vecb(:,4)=k2329*(Uh2-Phh+Un2.^2-2*U.*Un2+Phh34)-k83*Vhn;
        vecb(1,4)=vecb(1,4)+k2329*(-Chu1+Chu341);
        vecb(N-1,4)=vecb(N-1,4)+k2329*(-ChuN+Chu34N);
    
        U=U+phipm(k,J,vecb,10^(-10),1,1);

        t=t+k;        
             
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
    % calculated at the beginning of the next iteration
    n=2*n;   
end
