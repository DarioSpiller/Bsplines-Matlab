function [x,xdot,xddot]=Traiettoria_dot_off(A,N,T,dTdtau,ddTddtau,tau,dtaudt) 

%global A

% A= matrice dei coefficienti 
% N=numero di polinomi usati =D
% n=numeri di punti usati per interpolare la cirva


%tau=2*(t-t_i)/(t_f_ad-t_i)-1;
% dtaudt=2/(t_f-t_i);
% tau=linspace(-1,1,n);
% T(1,:)=ones(1,n);
% T(2,:)=tau;
% 
% dTdtau(1,:)=zeros(1,n);
% dTdtau(2,:)=ones(1,n);
% 
% ddTddtau(1,:)=zeros(1,n);
% ddTddtau(2,:)=zeros(1,n);
% ddTddtau(3,:)=4*ones(1,n);
% ddTddtau(4,:)=24*tau;
% 
% U(1,:)=ones(1,n);
% U(2,:)=2*tau;
% dTdtau2=dTdtau;


for ii=3:N
    
%     U(ii,:)=2*tau.*U(ii-1,:)-U(ii-2,:);  %polinomio di Chebyshev di II specie
    
    T(ii,:)=2*tau.*T(ii-1,:)-T(ii-2,:);  %traietoria
    
%       dTdtau(ii,:)=(ii-1)*U(ii-1,:);    %velocità
      dTdtau(ii,:)=2*T(2,:).*dTdtau(ii-1,:)+2*T(ii-1,:)-dTdtau(ii-2,:);
    if ii>4
        
   ddTddtau(ii,:)= (ii-1)*(2*dTdtau(ii-1,:)+ddTddtau(ii-2,:)/(ii-3));
    end

    
end

x=A*T;

xdot=A*dTdtau*dtaudt;

xddot=A*ddTddtau*dtaudt^2;%dtaudt^2;

end