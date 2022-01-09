close all;clear all;clc

% nodes' coordinates
N = [1,0,0; 2,0,40; 3,40,0; 4,40,40; 5,80,0; 6,80,40];

E = 1e7; % young's modulus
A = 1.5; % cross-sectional area

% element-node connectivity table
EN = [1,1,3; 2,1,4; 3,2,4; 4,3,4; 5,3,5; 6,5,4; 7,4,6; 8,5,6];

% length of the elements
for i=1:max(EN(:,1))
    L(i) = sqrt((N(EN(i,3),3)-N(EN(i,2),3)).^2+(N(EN(i,3),2)-N(EN(i,2),2)).^2);
end
L=L';
k=A*E./L;

% stiff matrix for each element
for i=1:max(EN(:,1))
    c(i) = (N(EN(i,3),2)-N(EN(i,2),2))./L(i);
    s(i) = (N(EN(i,3),3)-N(EN(i,2),3))./L(i);
    m = [c(i)^2,s(i)*c(i);s(i)*c(i),s(i)^2];
    K_e(:,:,i) = k(i).*[m,-m; -m,m];
end

% element global displacement location vector for each element
L_e = [2*EN(:,2)-1, 2*EN(:,2), 2*EN(:,3)-1, 2*EN(:,3)];

% global matrix
K = zeros(2*max(N(:,1)),2*max(N(:,1)));
for i=1:max(EN(:,1))
    K(L_e(i,:),L_e(i,:)) = K_e(:,:,i)+K(L_e(i,:),L_e(i,:));
end

% computes displacements and reaction forces
% nodes 1 and 2 are fixed, i.e U1=U2=U3=U4=0
U = zeros(2*max(N(:,1)),1);
F = zeros(2*max(N(:,1)),1);
U(1)=0; U(2)=0; U(3)=0; U(4)=0;
F(6)=-2000; F(9)=2000; F(11)=4000; F(12)=6000; %external forces

K_cc = K(1:4,1:4);
K_ca = K(1:4,5:12);
K_aa = K(5:12,5:12);

F_a = F(5:12);
U_c = U(1:4); % constrained displacements

U_a = K_aa^(-1)*(F_a-transpose(K_ca)*U_c); % active displacements
U(5:12) = U_a;

F_c = K_cc*U_c+K_ca*U_a; % reaction forces
F(1:4) = F_c;

% computes strains and stresses
for i=1:max(EN(:,1))
    eps(i) = (U(2*EN(i,3)-1)*c(i)+U(2*EN(i,3))*s(i)...
            -(U(2*EN(i,2)-1)*c(i)+U(2*EN(i,2))*s(i)))/L(i); % strains
    sig(i) = E*eps(i); % stresses
end
sig=sig'; eps=eps';

% displacements and forces
for i=1:max(N(:,1))
    U_x(i) = U(2*i-1);
    U_y(i) = U(2*i);
    F_x(i) = F(2*i-1);
    F_y(i) = F(2*i);
end
U_x=U_x'; U_y=U_y'; F_x=F_x'; F_y=F_y';

% results
A1 = [EN(:,1),eps,sig];
A2 = [N(:,1),U_x,U_y,F_x,F_y];
T1 = array2table(A1,...
    'VariableNames',{'Element','Strain','Stress (psi)'});
T2 = array2table(A2,...
    'VariableNames',{'Node','x-displacement (in)','y-displacement (in)','F_x (lb)','F_y (lb)'});
disp(T1); disp(T2);
