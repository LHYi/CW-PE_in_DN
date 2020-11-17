clear
clc
 
V0=1.02;

ZL01=0.01;
ZL12=0.03+1i*0.01;
ZL23=0.05+1i*0.01;
Zload1=1/(0.1-1i*0.03);
Sgen2=0.01;
Sload3=0.8-1i*0.01;

YL01=1/ZL01;
YL12=1/ZL12;
YL23=1/ZL23;
Yload1=1/Zload1;
Pgen2=real(Sgen2);
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
 
f=@(V1,V2,V3,theta1,theta2,theta3) [V0*V1*(G01*cos(theta1)+B01*sin(theta1))+V1*V1*G11+V1*V2*(G12*cos(theta1-theta2)+B12*sin(theta1-theta2));     V1*V2*(G12*cos(theta2-theta1)+B12*sin(theta2-theta1))+V2*V2*G22+V2*V3*(G23*cos(theta2-theta3)+B23*sin(theta2-theta3))-Pgen2;    V2*V3*(G23*cos(theta3-theta2)+B23*sin(theta3-theta2))+V3*V3*G33+Pload3;    V0*V1*(G01*sin(theta1)-B01*cos(theta1))-V1*V1*B11+V1*V2*(G12*sin(theta1-theta2)-B12*cos(theta1-theta2));    V1*V2*(G12*sin(theta2-theta1)-B12*cos(theta2-theta1))-V2*V2*B22+V2*V3*(G23*sin(theta2-theta3)-B23*cos(theta2-theta3));    V2*V3*(G23*sin(theta3-theta2)-B23*cos(theta3-theta2))-V3*V3*B33+Qload3];
 
J=@(V1,V2,V3,theta1,theta2,theta3) [V0*(G01*cos(theta1)+B01*sin(theta1))+2*V1*G11+V2*(G12*cos(theta1-theta2)+B12*sin(theta1-theta2)),    V1*(G12*cos(theta1-theta2)+B12*sin(theta1-theta2)),    0,    V0*V1*(-G01*sin(theta1)+B01*cos(theta1))+V1*V2*(-G12*sin(theta1-theta2)+B12*cos(theta1-theta2)),    V1*V2*(G12*sin(theta1-theta2)-B12*cos(theta1-theta2)),    0;    V2*(G12*cos(theta2-theta1)+B12*sin(theta2-theta1)),    V1*(G12*cos(theta2-theta1)+B12*sin(theta2-theta1))+2*V2*G22+V3*(G23*cos(theta2-theta3)+B23*sin(theta2-theta3)),    V2*(G23*cos(theta2-theta3)+B23*sin(theta2-theta3)),    V1*V2*(G12*sin(theta2-theta1)-B12*cos(theta2-theta1)),    V1*V2*(-G12*sin(theta2-theta1)+B12*cos(theta2-theta1))+V2*V3*(-G12*sin(theta2-theta3)+B23*cos(theta2-theta3)),    V2*V3*(G23*sin(theta2-theta3)-B23*cos(theta2-theta3));    0,    V3*(G23*cos(theta3-theta2)+B23*sin(theta3-theta2)),    V2*(G23*cos(theta3-theta2)+B23*sin(theta3-theta2))+2*V3*G33,    0,    V2*V3*(G23*sin(theta3-theta2)-B23*cos(theta3-theta2)),    V2*V3*(-G23*sin(theta3-theta2)+B23*cos(theta3-theta2));    V0*(G01*sin(theta1)-B01*cos(theta1))-2*V1*B11+V2*(G12*sin(theta1-theta2)-B12*cos(theta1-theta2)),    V1*(G12*sin(theta1-theta2)-B12*cos(theta1-theta2)),    0,    V0*V1*(G01*cos(theta1)+B01*sin(theta1))+V1*V2*(G12*cos(theta1-theta2)+B12*sin(theta1-theta2)),    V1*V2*(-G12*cos(theta1-theta2)-B12*sin(theta1-theta2)),    0;    V2*(G12*sin(theta2-theta1)-B12*cos(theta2-theta1)),    V1*(G12*sin(theta2-theta1)-B12*cos(theta2-theta1))-2*V2*B22+V3*(G23*sin(theta2-theta3)-B23*cos(theta2-theta3)),    V2*(G23*sin(theta2-theta3)-B23*cos(theta2-theta3)),    V1*V2*(-G12*cos(theta2-theta1)-B12*sin(theta2-theta1)),    V1*V2*(G12*cos(theta2-theta1)+B12*sin(theta2-theta1))+V2*V3*(G23*cos(theta2-theta3)+B23*sin(theta2-theta3)),    V2*V3*(-G23*cos(theta2-theta3)-B23*sin(theta2-theta3));    0,    V3*(G23*sin(theta3-theta2)-B23*cos(theta3-theta2)),    V2*(G23*sin(theta3-theta2)-B23*cos(theta3-theta2))-2*V3*B33,    0,    V2*V3*(-G23*cos(theta3-theta2)-B23*sin(theta3-theta2)),    V2*V3*(G23*cos(theta3-theta2)+B23*sin(theta3-theta2))];

fp=@(x) f(x(1),x(2),x(3),x(4),x(5),x(6));
Jp=@(x) J(x(1),x(2),x(3),x(4),x(5),x(6));
 
x=zeros(3,6);
x(1,:)=[1.02,1.02,1.02,0,0,0];
disp(norm(fp(x(1,:))));
for n=2:length(x)
    x(n,:)=x(n-1,:)-((Jp(x(n-1,:))^(-1))*fp(x(n-1,:)))';
    disp(norm(fp(x(n,:))));
end
disp('Solution:');
disp(x(n,:));