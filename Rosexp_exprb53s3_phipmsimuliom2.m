% Program that implements exponential Rosenbrock method given by tableau
%
% 0      |   
% 1/2    |  a_21
% 9/10   |  a_31   a_32    
% -----------------------------------
%        |  b_1     b_2      b_3        
% with 
% a_21=1/2 phi_{1,2}
% a_31=9/10 phi_{1,3}-27/25 phi_{3,2}-729/125 phi_{3,3}
% a_32=27/25 phi_{3,2}+729/125 phi_{3,3}
% b_1=phi_1-1208/81 phi3+1120/27 phi4  
% b_2=18 phi3-60 phi_4
% b_3=-250/81 phi_3+500/27 phi_4
%
% with the standard method of lines 
% u_t=u_xx+f(t,u), x\in [0,1]
% u(t,0)=g0(t), u(t,1)=g1(t) 
% with f(t,u)=u^2+h(t,u) and with exact solution u(t,x)=cos(t+x)
% For the spatial discretization we have used the standard second-order
% difference squeme. 
% This methods appears in 
% V. T. Luan and A. Ostermann, Exponential Rosenbrock methods of order 
% five-construction, analysis and numerical comparisons, J. Comput. Appl. 
% Math. 255 (2014) 417-431
% Here, we use subrutine phipm_simul_iom in
% V. T. Luan, J. A. Pudykiewicz and D. R. Reynolds, Further development of
% efficient and accurate time integration schemes for meteorological
% models, J. Comp. Phys. 376 (2019) 817-837

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
    k12=k/2;
    k910=9*k/10;
    k2=2*k;
    k569=56/(9*k);
    k22=k^2;
    k2281=81*k22;
    k22403=40/(k22*3);
    k33=k^3;
    k3327=k33*27;

    tie(1)=k12;
    tie(2)=k910;
    tie(3)=k;

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
    Phh910=zeros(N-1,1);
    Vhn=zeros(N-1,1);
    Un2=zeros(N-1,1);
    Un3=zeros(N-1,1);
    Gn1=zeros(N-1,1);
    Gn2=zeros(N-1,1);
    vecb=zeros(N-1,5);
    vecb2=zeros(N-1,3);
    vecb31=zeros(N-1,4);
    vecb32=zeros(N-1,4);

    % CPU time calculus starts 
    tstart=tic;

    % For the local error r=1. For the global one, r=n
    for r=1:n
        % J is the Jacobian
        J=A+2*diag(U);
        J=sparse(J);        
        t12=t+k12;
        t910=t+k910;
        
        % Boundary values that are needed. They are 0 except for ii=1 and ii=N.
        % They are Dirichlet. uu calculates the boundary values of u and uut
        % calculates the boundary values of u_t            
        Chu1=hdiv*uu(t,0);
        ChuN=hdiv*uu(t,1);
    
        Chu121=hdiv*uu(t12,0);
        Chu12N=hdiv*uu(t12,1);  
    
        Chu9101=hdiv*uu(t910,0);
        Chu910N=hdiv*uu(t910,1);      
    
        Chut1=hdiv*uut(t,0);
        ChutN=hdiv*uut(t,1);
        
        % funh is funcion h(t,x) and function funht is h_t(t,x)           
        for ii=1:N-1
            Phh(ii)=funh(t,x(ii)); 
            Phh12(ii)=funh(t12,x(ii));
            Phh910(ii)=funh(t910,x(ii));
            Vhn(ii)=funht(t,x(ii));        
        end
    
        Vhn(1)=Vhn(1)+Chut1;
        Vhn(N-1)=Vhn(N-1)+ChutN;  
        
        % vecb and vecb2 contain the vectors that are multiplied by the
        % exponential funcions. The first column is zero, the second one is
        % multiplied by phi_1(c_ikA), the third one by phi_2(c_ikA)...

        % Stages
        Au=A*U;
        Uh2=U.^2;
        vecb2(:,2)=Au+Uh2+Phh;
        vecb2(1,2)=vecb2(1,2)+Chu1;
        vecb2(N-1,2)=vecb2(N-1,2)+ChuN;
        
        vecb2(:,3)=Vhn; 
    
        Uaux=phipm_simul_iom(tie,J,vecb2,10^(-10),1,3);

        Un2=U+Uaux(:,1);

        % Second stage. The part that is multiplied by phi3(9/10) and
        % phi3(1/2) only diffrs in a constant 
        Uh2Ph=-U.^2+Phh;
        Un2Ph=Un2.^2+Phh12-2*U.*Un2;   
        vecb32(:,4)=(-Uh2Ph+Un2Ph)/k22-Vhn/k2;
        vecb32(1,4)=vecb32(1,4)+(Chu121-Chu1)/k22;    
        vecb32(N-1,4)=vecb32(N-1,4)+(Chu12N-ChuN)/k22;  

        U3aux=phipm_simul_iom(tie(1,1:2),J,vecb32,10^(-10),1,2);
    
        Un3=U+Uaux(:,2)+8*U3aux(:,1)+216*U3aux(:,2)/25; 
 
        % Aproximation to the exact solution at time t_{n+1}  with 
        % horizontal/vertical implementation             
        Un3Ph=Un3.^2+Phh910-2*U.*Un3;
        vecb(:,4)=(1458*Un2Ph-1208*Uh2Ph-250*Un3Ph)/k2281-k569*Vhn;
        vecb(1,4)=vecb(1,4)+(1458*Chu121-1208*Chu1-250*Chu9101)/k2281;
        vecb(N-1,4)=vecb(N-1,4)+(1458*Chu12N-1208*ChuN-250*Chu910N)/k2281;  
    
        vecb(:,5)=(-1620*Un2Ph+1120*Uh2Ph+500*Un3Ph)/k3327+k22403*Vhn;
        vecb(1,5)=vecb(1,5)+(-1620*Chu121+1120*Chu1+500*Chu9101)/k3327;
        vecb(N-1,5)=vecb(N-1,5)+(-1620*Chu12N+1120*ChuN+500*Chu910N)/k3327;
    
        U=U+Uaux(:,3)+phipm(k,J,vecb,10^(-10),1,1); 
    
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
    % calculated at the beginning os the next iteration
    n=2*n;

end
