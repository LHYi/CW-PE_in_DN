clear
clc
 
V0=1.02;

ZL01=0.01;
ZL12=0.03+1i*0.01;
ZL23=0.05+1i*0.01;
Zload1=1/(0.1-1i*0.03);
Sload3=0.8-1i*0.01;

YL01=1/ZL01;
YL12=1/ZL12;
Yload1=1/Zload1;
Pload2=0.79;
Qload2=0;
 
Y=[YL01,-YL01,0;
-YL01,YL01+YL12+Yload1,-YL12;
0,-YL12,YL12;];

G01=real(Y(1,2));
B01=imag(Y(1,2));
 
G11=real(Y(2,2));
B11=imag(Y(2,2));
 
G12=real(Y(2,3));
B12=imag(Y(2,3));
 
G22=real(Y(3,3));
B22=imag(Y(3,3));
 
k=100;
V1 = zeros(1,k);
V2 = zeros(1,k);
baseline = ones(1,k);
%Q = zeros(1,k);
P = zeros(1,k);
Loss01=zeros(1,k);
Loss12=zeros(1,k);
Loss=zeros(1,k);
for m = 1:k
    Pload2 = Pload2+0.01;
    f=@(V1,V2,theta1,theta2) [V0*V1*(G01*cos(theta1)+B01*sin(theta1))+V1*V1*G11+V1*V2*(G12*cos(theta1-theta2)+B12*sin(theta1-theta2));     V1*V2*(G12*cos(theta2-theta1)+B12*sin(theta2-theta1))+V2*V2*G22+Pload2;    V0*V1*(G01*sin(theta1)-B01*cos(theta1))-V1*V1*B11+V1*V2*(G12*sin(theta1-theta2)-B12*cos(theta1-theta2));    V1*V2*(G12*sin(theta2-theta1)-B12*cos(theta2-theta1))-V2*V2*B22-Qload2;];

    J=@(V1,V2,theta1,theta2) [V0*(G01*cos(theta1)+B01*sin(theta1))+2*V1*G11+V2*(G12*cos(theta1-theta2)+B12*sin(theta1-theta2)),    V1*(G12*cos(theta1-theta2)+B12*sin(theta1-theta2)),    V0*V1*(-G01*sin(theta1)+B01*cos(theta1))+V1*V2*(-G12*sin(theta1-theta2)+B12*cos(theta1-theta2)),    V1*V2*(G12*sin(theta1-theta2)-B12*cos(theta1-theta2));    V2*(G12*cos(theta2-theta1)+B12*sin(theta2-theta1)),    V1*(G12*cos(theta2-theta1)+B12*sin(theta2-theta1))+2*V2*G22,    V1*V2*(G12*sin(theta2-theta1)-B12*cos(theta2-theta1)),    V1*V2*(-G12*sin(theta2-theta1)+B12*cos(theta2-theta1));    V0*(G01*sin(theta1)-B01*cos(theta1))-2*V1*B11+V2*(G12*sin(theta1-theta2)-B12*cos(theta1-theta2)),    V1*(G12*sin(theta1-theta2)-B12*cos(theta1-theta2)),    V0*V1*(G01*cos(theta1)+B01*sin(theta1))+V1*V2*(G12*cos(theta1-theta2)+B12*sin(theta1-theta2)),    V1*V2*(-G12*cos(theta1-theta2)-B12*sin(theta1-theta2));    V2*(G12*sin(theta2-theta1)-B12*cos(theta2-theta1)),    V1*(G12*sin(theta2-theta1)-B12*cos(theta2-theta1))-2*V2*B22,    V1*V2*(-G12*cos(theta2-theta1)-B12*sin(theta2-theta1)),    V1*V2*(G12*cos(theta2-theta1)+B12*sin(theta2-theta1));];

    fp=@(x) f(x(1),x(2),x(3),x(4));
    Jp=@(x) J(x(1),x(2),x(3),x(4));

    x=zeros(3,4);
    x(1,:)=[1.02,1.02,0,0];
    %disp(norm(fp(x(1,:))));
    for n=2:length(x)
        x(n,:)=x(n-1,:)-((Jp(x(n-1,:))^(-1))*fp(x(n-1,:)))';
        %disp(norm(fp(x(n,:))));
    end
    %disp('Solution:');
    %disp(x(n,:));
    V1(m)=x(n,1);
    V2(m)=x(n,2);
    P(m)=Pload2;
    Loss12(m)=(Pload2*Pload2+Qload2*Qload2)/(V2(m)*V2(m))*ZL12;
    S1=V1(m)*V1(m)*Yload1+Loss12(m)+Pload2+1i*Qload2;
    P1=real(S1);
    Q1=imag(S1);
    Loss01(m)=(P1*P1+Q1*Q1)/(V1(m)*V1(m))*ZL01;
    Loss(m)=Loss01(m)+Loss12(m);
end
figure;
plot(P,V1,'b',P,V2,'r');
title('Node Voltage changes as P3 changes');
xlabel('Active power transferred by MFC (p.u.)');
ylabel('Voltage at the nodes (p.u.)');
legend('V1','V2','V=1');
% figure;
% plot(Q,Loss,Q,Loss01,Q,Loss12);
% hold on;
% plot(Q(find(Loss==min(Loss))),min(Loss),'r*');
