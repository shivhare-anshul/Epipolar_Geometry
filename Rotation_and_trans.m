function [R1,t]= Rotation_and_trans(E,x)
 
x(1,3) = 1; x=x';
W1 = [0 -1 0;1 0 0;0 0 1];
W2 = W1';
[U,~,V] = svd(E);
%S=S/S(2,2);
R1 = (-1)*U * W1 * V';
R2 = (-1)*U * W2 * V';
R3 = U * W1 * V';
R4 = U * W2 * V';
t = U(:,3);
t = t/norm(t);

P1 = zeros(3,4); P2 = zeros(3,4); P3 = zeros(3,4); P4 = zeros(3,4);

P1(1:3,1:3) = R1;  P1(1:3,4)= t ;
P2(1:3,1:3) = R2;  P2(1:3,4)= t ;
P3(1:3,1:3) = R3;  P3(1:3,4)= t ;
P4(1:3,1:3) = R4;  P4(1:3,4)= t ;

a1= pinv(P1)*x; a1 = a1/a1(4,1)
a2 = pinv(P2)*x;  a2 = a2/a2(4,1)
a3 = pinv(P3)*x;  a3 = a3/a3(4,1)
a4 = pinv(P4)*x;  a4 = a4/a4(4,1)

det(R1)
det(R2)
det(R3)
det(R4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Z=[0 1 0;-1 0 0;0 0 0];
% S = U*Z*U';
% t=null(S);
% t = t/norm(t);
end
