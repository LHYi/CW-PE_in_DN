clear
clc
 
V0=1.02;

ZL01=0.01;
ZL12=0.03+1i*0.01;
Zload1=1/(0.1-1i*0.03);
Sload2=0.8347;
Sgen2=0.01;
Pgen2=real(Sgen2);

YL01=1/ZL01;
YL12=1/ZL12;
Yload1=1/Zload1;
Pload2=real(Sload2);
 
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
 
f=@(V1,V2,theta1,theta2) [V0*V1*(G01*cos(theta1)+B01*sin(theta1))+V1*V1*G11+V1*V2*(G12*cos(theta1-theta2)+B12*sin(theta1-theta2));     V1*V2*(G12*cos(theta2-theta1)+B12*sin(theta2-theta1))+V2*V2*G22+Pload2+Pgen2;    V0*V1*(G01*sin(theta1)-B01*cos(theta1))-V1*V1*B11+V1*V2*(G12*sin(theta1-theta2)-B12*cos(theta1-theta2));    V1*V2*(G12*sin(theta2-theta1)-B12*cos(theta2-theta1))-V2*V2*B22;];
 
J=@(V1,V2,theta1,theta2) [V0*(G01*cos(theta1)+B01*sin(theta1))+2*V1*G11+V2*(G12*cos(theta1-theta2)+B12*sin(theta1-theta2)),    V1*(G12*cos(theta1-theta2)+B12*sin(theta1-theta2)),    V0*V1*(-G01*sin(theta1)+B01*cos(theta1))+V1*V2*(-G12*sin(theta1-theta2)+B12*cos(theta1-theta2)),    V1*V2*(G12*sin(theta1-theta2)-B12*cos(theta1-theta2));    V2*(G12*cos(theta2-theta1)+B12*sin(theta2-theta1)),    V1*(G12*cos(theta2-theta1)+B12*sin(theta2-theta1))+2*V2*G22,    V1*V2*(G12*sin(theta2-theta1)-B12*cos(theta2-theta1)),    V1*V2*(-G12*sin(theta2-theta1)+B12*cos(theta2-theta1));    V0*(G01*sin(theta1)-B01*cos(theta1))-2*V1*B11+V2*(G12*sin(theta1-theta2)-B12*cos(theta1-theta2)),    V1*(G12*sin(theta1-theta2)-B12*cos(theta1-theta2)),    V0*V1*(G01*cos(theta1)+B01*sin(theta1))+V1*V2*(G12*cos(theta1-theta2)+B12*sin(theta1-theta2)),    V1*V2*(-G12*cos(theta1-theta2)-B12*sin(theta1-theta2));    V2*(G12*sin(theta2-theta1)-B12*cos(theta2-theta1)),    V1*(G12*sin(theta2-theta1)-B12*cos(theta2-theta1))-2*V2*B22,    V1*V2*(-G12*cos(theta2-theta1)-B12*sin(theta2-theta1)),    V1*V2*(G12*cos(theta2-theta1)+B12*sin(theta2-theta1));];

fp=@(x) f(x(1),x(2),x(3),x(4));
Jp=@(x) J(x(1),x(2),x(3),x(4));
 
x=zeros(3,4);
x(1,:)=[1.02,1.02,0,0];
disp(norm(fp(x(1,:))));
for n=2:length(x)
    x(n,:)=x(n-1,:)-((Jp(x(n-1,:))^(-1))*fp(x(n-1,:)))';
    disp(norm(fp(x(n,:))));
end
disp('Solution:');
disp(x(n,:));