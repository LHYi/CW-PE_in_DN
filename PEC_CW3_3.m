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
YL23=1/ZL23;
Yload1=1/Zload1;
Pload3=real(Sload3);
Qload3=imag(Sload3);
 
Y=[YL01,-YL01,0,0;
-YL01,YL01+YL12+Yload1,-YL12,0;
0,-YL12,YL12+YL23,-YL23;
0,0,-YL23,YL23];

G01=real(Y(1,2));
B01=imag(Y(1,2));
 
G11=real(Y(2,2));
B11=imag(Y(2,2));
 
G12=real(Y(2,3));
B12=imag(Y(2,3));
 
G22=real(Y(3,3));
B22=imag(Y(3,3));

G23=real(Y(3,4));
B23=imag(Y(3,4));

G33=real(Y(4,4));
B33=imag(Y(4,4));

k=100;
V1 = zeros(1,k);
V2 = zeros(1,k);
V2b = zeros(1,k);
V3 = zeros(1,k);
%Q = zeros(1,k);
P3 = zeros(1,k);
Loss01=zeros(1,k);
Loss12=zeros(1,k);
Loss=zeros(1,k);
Qload2=0;

for m = 1:k
    V3=0.96;
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
    f1=@(V2b,theta3) [V2b*V3*(G23*cos(theta3)+B23*sin(theta3))+V3*V3*G33+Pload3; V2b*V3*(G23*sin(theta3)-B23*cos(theta3))-V3*V3*B33+Qload3;];
    J1=@(V2b,theta3) [V3*(G23*cos(theta3)+B23*sin(theta3)), V2b*V3*(-G23*sin(theta3)+B23*cos(theta3)); V3*(G23*sin(theta3)-B23*cos(theta3)), V2b*V3*(G23*cos(theta3)+B23*sin(theta3));];

    fp1=@(x1) f1(x1(1),x1(2));
    Jp1=@(x1) J1(x1(1),x1(2));
 
    x1=zeros(5,2);
    x1(1,:)=[1.02,0];
    %disp(norm(fp(x1(1,:))));
    for n=2:length(x1)
        x1(n,:)=x1(n-1,:)-((Jp1(x1(n-1,:))^(-1))*fp1(x1(n-1,:)))';
        %disp(norm(fp1(x1(n,:))));
    end
    
    V2b(m) = x1(n,1);
    theta3 = x1(n,2);
    P3(m)=Pload3;
    Pload2 = V2b(m)*V3*(G23*cos(-theta3)+B23*sin(-theta3))+V2b(m)*V2b(m)*G33;
    
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
    
    f2=@(V1,V2,theta1,theta2) [V0*V1*(G01*cos(theta1)+B01*sin(theta1))+V1*V1*G11+V1*V2*(G12*cos(theta1-theta2)+B12*sin(theta1-theta2));     V1*V2*(G12*cos(theta2-theta1)+B12*sin(theta2-theta1))+V2*V2*G22+Pload2;    V0*V1*(G01*sin(theta1)-B01*cos(theta1))-V1*V1*B11+V1*V2*(G12*sin(theta1-theta2)-B12*cos(theta1-theta2));    V1*V2*(G12*sin(theta2-theta1)-B12*cos(theta2-theta1))-V2*V2*B22-Qload2;];
    J2=@(V1,V2,theta1,theta2) [V0*(G01*cos(theta1)+B01*sin(theta1))+2*V1*G11+V2*(G12*cos(theta1-theta2)+B12*sin(theta1-theta2)),    V1*(G12*cos(theta1-theta2)+B12*sin(theta1-theta2)),    V0*V1*(-G01*sin(theta1)+B01*cos(theta1))+V1*V2*(-G12*sin(theta1-theta2)+B12*cos(theta1-theta2)),    V1*V2*(G12*sin(theta1-theta2)-B12*cos(theta1-theta2));    V2*(G12*cos(theta2-theta1)+B12*sin(theta2-theta1)),    V1*(G12*cos(theta2-theta1)+B12*sin(theta2-theta1))+2*V2*G22,    V1*V2*(G12*sin(theta2-theta1)-B12*cos(theta2-theta1)),    V1*V2*(-G12*sin(theta2-theta1)+B12*cos(theta2-theta1));    V0*(G01*sin(theta1)-B01*cos(theta1))-2*V1*B11+V2*(G12*sin(theta1-theta2)-B12*cos(theta1-theta2)),    V1*(G12*sin(theta1-theta2)-B12*cos(theta1-theta2)),    V0*V1*(G01*cos(theta1)+B01*sin(theta1))+V1*V2*(G12*cos(theta1-theta2)+B12*sin(theta1-theta2)),    V1*V2*(-G12*cos(theta1-theta2)-B12*sin(theta1-theta2));    V2*(G12*sin(theta2-theta1)-B12*cos(theta2-theta1)),    V1*(G12*sin(theta2-theta1)-B12*cos(theta2-theta1))-2*V2*B22,    V1*V2*(-G12*cos(theta2-theta1)-B12*sin(theta2-theta1)),    V1*V2*(G12*cos(theta2-theta1)+B12*sin(theta2-theta1));];

    fp2=@(x2) f2(x2(1),x2(2),x2(3),x2(4));
    Jp2=@(x2) J2(x2(1),x2(2),x2(3),x2(4));
    
    x2=zeros(5,4);
    x2(1,:)=[1.02,1.02,0,0];
    %disp(norm(fp(x(1,:))));
    for n=2:length(x2)
        x2(n,:)=x2(n-1,:)-((Jp2(x2(n-1,:))^(-1))*fp2(x2(n-1,:)))';
        %disp(norm(fp(x(n,:))));
    end
    %disp('Solution:');
    %disp(x(n,:));
    V1(m)=x2(n,1);
    V2(m)=x2(n,2);
    Loss12(m)=(Pload2*Pload2+Qload2*Qload2)/(V2(m)*V2(m))*ZL12;
    S1=V1(m)*V1(m)*Yload1+Loss12(m)+Pload2+1i*Qload2;
    P1=real(S1);
    Q1=imag(S1);
    Loss01(m)=(P1*P1+Q1*Q1)/(V1(m)*V1(m))*ZL01;
    Loss(m)=Loss01(m)+Loss12(m);
    Pload3 = Pload3+0.01;
end

figure;
plot(P3,V1,'b',P3,V2,'r',P3,V2b,'m');
hold on;

V1 = zeros(1,k);
V2 = zeros(1,k);
V2b = zeros(1,k);
V3 = zeros(1,k);
%Q = zeros(1,k);
P3 = zeros(1,k);
Loss01=zeros(1,k);
Loss12=zeros(1,k);
Loss=zeros(1,k);
Qload2=0.5;
Pload3=real(Sload3);
for m = 1:k
    V3=0.96;
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
    f1=@(V2b,theta3) [V2b*V3*(G23*cos(theta3)+B23*sin(theta3))+V3*V3*G33+Pload3; V2b*V3*(G23*sin(theta3)-B23*cos(theta3))-V3*V3*B33+Qload3;];
    J1=@(V2b,theta3) [V3*(G23*cos(theta3)+B23*sin(theta3)), V2b*V3*(-G23*sin(theta3)+B23*cos(theta3)); V3*(G23*sin(theta3)-B23*cos(theta3)), V2b*V3*(G23*cos(theta3)+B23*sin(theta3));];

    fp1=@(x1) f1(x1(1),x1(2));
    Jp1=@(x1) J1(x1(1),x1(2));
 
    x1=zeros(5,2);
    x1(1,:)=[1.02,0];
    %disp(norm(fp(x1(1,:))));
    for n=2:length(x1)
        x1(n,:)=x1(n-1,:)-((Jp1(x1(n-1,:))^(-1))*fp1(x1(n-1,:)))';
        %disp(norm(fp1(x1(n,:))));
    end
    
    V2b(m) = x1(n,1);
    theta3 = x1(n,2);
    P3(m)=Pload3;
    Pload2 = V2b(m)*V3*(G23*cos(-theta3)+B23*sin(-theta3))+V2b(m)*V2b(m)*G33;
    
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
    
    f2=@(V1,V2,theta1,theta2) [V0*V1*(G01*cos(theta1)+B01*sin(theta1))+V1*V1*G11+V1*V2*(G12*cos(theta1-theta2)+B12*sin(theta1-theta2));     V1*V2*(G12*cos(theta2-theta1)+B12*sin(theta2-theta1))+V2*V2*G22+Pload2;    V0*V1*(G01*sin(theta1)-B01*cos(theta1))-V1*V1*B11+V1*V2*(G12*sin(theta1-theta2)-B12*cos(theta1-theta2));    V1*V2*(G12*sin(theta2-theta1)-B12*cos(theta2-theta1))-V2*V2*B22-Qload2;];
    J2=@(V1,V2,theta1,theta2) [V0*(G01*cos(theta1)+B01*sin(theta1))+2*V1*G11+V2*(G12*cos(theta1-theta2)+B12*sin(theta1-theta2)),    V1*(G12*cos(theta1-theta2)+B12*sin(theta1-theta2)),    V0*V1*(-G01*sin(theta1)+B01*cos(theta1))+V1*V2*(-G12*sin(theta1-theta2)+B12*cos(theta1-theta2)),    V1*V2*(G12*sin(theta1-theta2)-B12*cos(theta1-theta2));    V2*(G12*cos(theta2-theta1)+B12*sin(theta2-theta1)),    V1*(G12*cos(theta2-theta1)+B12*sin(theta2-theta1))+2*V2*G22,    V1*V2*(G12*sin(theta2-theta1)-B12*cos(theta2-theta1)),    V1*V2*(-G12*sin(theta2-theta1)+B12*cos(theta2-theta1));    V0*(G01*sin(theta1)-B01*cos(theta1))-2*V1*B11+V2*(G12*sin(theta1-theta2)-B12*cos(theta1-theta2)),    V1*(G12*sin(theta1-theta2)-B12*cos(theta1-theta2)),    V0*V1*(G01*cos(theta1)+B01*sin(theta1))+V1*V2*(G12*cos(theta1-theta2)+B12*sin(theta1-theta2)),    V1*V2*(-G12*cos(theta1-theta2)-B12*sin(theta1-theta2));    V2*(G12*sin(theta2-theta1)-B12*cos(theta2-theta1)),    V1*(G12*sin(theta2-theta1)-B12*cos(theta2-theta1))-2*V2*B22,    V1*V2*(-G12*cos(theta2-theta1)-B12*sin(theta2-theta1)),    V1*V2*(G12*cos(theta2-theta1)+B12*sin(theta2-theta1));];

    fp2=@(x2) f2(x2(1),x2(2),x2(3),x2(4));
    Jp2=@(x2) J2(x2(1),x2(2),x2(3),x2(4));
    
    x2=zeros(5,4);
    x2(1,:)=[1.02,1.02,0,0];
    %disp(norm(fp(x(1,:))));
    for n=2:length(x2)
        x2(n,:)=x2(n-1,:)-((Jp2(x2(n-1,:))^(-1))*fp2(x2(n-1,:)))';
        %disp(norm(fp(x(n,:))));
    end
    %disp('Solution:');
    %disp(x(n,:));
    V1(m)=x2(n,1);
    V2(m)=x2(n,2);
    Loss12(m)=(Pload2*Pload2+Qload2*Qload2)/(V2(m)*V2(m))*ZL12;
    S1=V1(m)*V1(m)*Yload1+Loss12(m)+Pload2+1i*Qload2;
    P1=real(S1);
    Q1=imag(S1);
    Loss01(m)=(P1*P1+Q1*Q1)/(V1(m)*V1(m))*ZL01;
    Loss(m)=Loss01(m)+Loss12(m);
    Pload3 = Pload3+0.01;
end
baseline=ones(1,k)*0.96;
plot(P3,V1,'--b',P3,V2,'--r',P3,V2b,'--m',P3,baseline,'k');
title('Node Voltage changes as Pload3 changes when QMFC=0 and QMFC=0.5');
xlabel('Pload3 (p.u.)');
ylabel('Voltage at the nodes (p.u.)');
legend('V1','V2','V2''','V1','V2','V2''');
lgd = legend;
lgd.NumColumns = 2;