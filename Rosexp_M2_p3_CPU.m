% Program that implements exponential Rosenbrock method given by tableau
%
% 0    |   
% 1    |  phi_1
% -------------------------------------------------------------------------
%      |  7/6 phi_1-phi2    phi2-1/6 phi_1    
%
% avoiding the order reduction with p=3 for problem 
% u_t=u_xx+f(t,u), x\in [0,1]
% u(t,0)=g0(t), u(t,1)=g1(t) 
% with f(t,u)=u^2+h(t,u) and with exact solution u(t,x)=cos(t+x)
% For the spatial discretization we have used the standard second-order
% difference squeme. 

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
    k22=k^2;
    k2212=k22/12;
    k2m=k/2;
    
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
    Phh1=zeros(N-1,1);
    Vhn=zeros(N-1,1);
    Un2=zeros(N-1,1);
    vecb=zeros(N-1,5);
    vecb2=zeros(N-1,4);

    % CPU time calculus starts 
    tstart=tic;
    
    % funh is funcion h(t,x)
    for ii=1:N-1
        Phh(ii)=funh(t,x(ii));         
    end
    
    % For the local error r=1. For the global one, r=n    
    for r=1:n
        % J is the Jacobian 
        J=A+2*diag(U);
        J=sparse(J);
        
        tk=t+k;

        % Boundary values that are needed. They are 0 except for ii=1 and ii=N.
        % They are Dirichlet. uu calculates the boundary values of u, uut
        % calculates the boundary values of u_t, uutt calculates the
        % boundaru values of u_tt, funhtt calculates function h_tt(t,x)
        Chu1=hdiv*uu(t,0);
        ChuN=hdiv*uu(t,1);
      
        aux1=uut(t,0);
        aux2=uut(t,1);
        Chut1=hdiv*aux1;
        ChutN=hdiv*aux2;

        Chutt1=hdiv*uutt(t,0);
        ChuttN=hdiv*uutt(t,1);
    
        % Boundary value of 2u_t^2+h_tt
        Chsum1=2*aux1*Chut1+hdiv*funhtt(t,0);
        ChsumN=2*aux2*ChutN+hdiv*funhtt(t,1);
        
        % funh is funcion h(t,x) and function funht is h_t(t,x)    
        for ii=1:N-1
            Phh1(ii)=funh(tk,x(ii));  
            Vhn(ii)=funht(t,x(ii));        
        end
    
        Uh2=U.^2;
        Au=A*U;
        
        % vecb2 and vecb contain the vectors that are multiplied by the
        % exponential funcions. The first column is zero, the second one is
        % multiplied by phi_1(c_ikA), the third one by phi_2(c_ikA)...
    
        % Stage
        vecb2(:,2)=Au+Uh2+Phh;
        vecb2(1,2)=vecb2(1,2)+Chu1;
        vecb2(N-1,2)=vecb2(N-1,2)+ChuN;

        vecb2(:,3)=Vhn;
        vecb2(1,3)=vecb2(1,3)+Chut1;
        vecb2(N-1,3)=vecb2(N-1,3)+ChutN;
    
        Un2=U+phipm(k,J,vecb2,10^(-10),1,1); 
        
        % Aproximation to the exact solution at time t_{n+1}  
        Uh2Ph=Un2.^2-2*U.*Un2+Phh1;  
        
        vecb(:,2)=Au+(5*Uh2+7*Phh-Uh2Ph+k*Vhn)/6;
        vecb(1,2)=vecb(1,2)+Chu1;
        vecb(N-1,2)=vecb(N-1,2)+ChuN;
    
        vecb(:,3)=(Uh2-Phh+Uh2Ph)/k;
        vecb(1,3)=vecb(1,3)+Chut1-k2212*Chsum1;
        vecb(N-1,3)=vecb(N-1,3)+ChutN-k2212*ChsumN;
       
        vecb(1,4)=Chutt1+k2m*Chsum1;            
        vecb(N-1,4)=ChuttN+k2m*ChsumN;   
        vecb(:,4)=sparse(vecb(:,4));
   
        % uuttt calculates the boundary values of u_ttt
        vecb(1,5)=hdiv*uuttt(t,0)-Chsum1;
        vecb(N-1,5)=hdiv*uuttt(t,1)-ChsumN;    
        vecb(:,5)=sparse(vecb(:,5));
   
        U=U+phipm(k,J,vecb,10^(-10),1,1);
    
        t=tk;        
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
    % calculated at the beginning of the next iteration
    n=2*n;

end
