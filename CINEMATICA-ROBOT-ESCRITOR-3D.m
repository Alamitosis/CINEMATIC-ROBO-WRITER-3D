%Programa de un robot manipuñador de GDL
%Trayectoria Recta parametrizada por polinomio de 5to grado
%Denavit-Hartenberg
%Dr. Waldemar Pérez Bailón
clc
clear all
  %%%%%%%%%%circulo de la a
l1=20;
l2=20;
cx=-10;
cy=5;
cz=0;
r=5;
largo=2*pi*r;
tf=20;
qf=largo;
 
t=0:0.1:tf;
q=qf*(10*(t.^3)/(tf^3)-15*(t.^4)/(tf^4)+6*(t.^5)/(tf^5));
 
for a=1:1:length(t)
    px(a)=r*sind(360*q(a)/qf)+cx;
    py(a)=r*cosd(360*q(a)/qf)+cy;
    pz(a)=0*r*cosd(360*q(a)/qf)+cz;
    
    theta1=atan2(py(a),px(a))*180/pi;
    C1=cosd(theta1);
    S1=sind(theta1);
    
    C3=((px(a)^2)+(py(a)^2)+(pz(a)^2)-(l2^2)-(l1^2))/(2*l1*l2);
    S3=(1-C3^2)^(1/2);
    theta3=atan2(S3,C3)*180/pi;
    
    C2=((l1+l2*C3)*(px(a)*C1+py(a)*S1)-l2*pz(a)*S3)/(pz(a)^2+(px(a)*C1+py(a)*S1)^2);
    S2=(1-C2^2)^(1/2);
    theta2=atan2(S2,C2)*180/pi;
 
    T01=[cosd(theta1)  ,-sind(theta1) ,0  ,0;
         sind(theta1)  ,cosd(theta1)  ,0  ,0;
         0             ,0             ,1  ,0;
         0             ,0             ,0  ,1];
 
    T12=[cosd(theta2)  ,-sind(theta2) ,0  ,0;
         0             ,0             ,-1 ,0;
         sind(theta2)  ,cosd(theta2)  ,0  ,0;
         0             ,0             ,0  ,1];
 
    T23=[cosd(theta3)  ,-sind(theta3) ,0  ,l1;
         -sind(theta3) ,-cosd(theta3) ,0  ,0;
         0             ,0             ,-1 ,0;
         0             ,0             ,0  ,1];
 
    T34=[1   ,0   ,0  ,l2;
         0   ,1   ,0  ,0;
         0   ,0   ,1  ,0;
         0   ,0   ,0  ,1];
 
    T02=T01*T12;
    T03=T02*T23;
    T04=T03*T34;
 
    Lx1=[T02(1,4) T03(1,4)];
    Ly1=[T02(2,4) T03(2,4)];
    Lz1=[T02(3,4) T03(3,4)];
    Lx2=[T03(1,4) T04(1,4)];
    Ly2=[T03(2,4) T04(2,4)];
    Lz2=[T03(3,4) T04(3,4)];
    
    plot3(Lx1,Ly1,Lz1,'r','LineWidth',1)
    grid on
    hold on
    axis([-20 40 -20 20 -0 20])
    plot3(Lx2,Ly2,Lz2,'b','LineWidth',1)
    hold on
    plot3(px,py,pz,'k','LineWidth',2)
    hold off
    pause(0.1)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%LINEA


%Programa de un robot manipuñador de GDL
%Trayectoria Circular parametrizada por polinomio de 5to grado
%Denavit-Hartenberg
%Dr. Waldemar Pérez Bailón

%%clc
%%clear all

l1=20;
l2=20;
 
p1=[-5 10 0];
p2=[-5 0 0];
 
x=p2(1)-p1(1);
y=p2(2)-p1(2);
z=p2(3)-p1(3);
 
largo=(x^2+y^2+z^2)^(1/2);
tf=50;
qf=largo;
t=0:0.1:tf;
q=qf*(10*(t.^3)/(tf^3)-15*(t.^4)/(tf^4)+6*(t.^5)/(tf^5));
 
 
for a=1:1:length(t)
    
    pX(a)=p1(1)+x*q(a)/largo;
    pY(a)=p1(2)+y*q(a)/largo;
    pZ(a)=p1(3)+z*q(a)/largo;
    
    theta1=atan2(pY(a),pX(a))*180/pi;
    C1=cosd(theta1);
    S1=sind(theta1);
    
    C3=((pX(a)^2)+(pY(a)^2)+(pZ(a)^2)-(l2^2)-(l1^2))/(2*l1*l2);
    S3=(1-C3^2)^(1/2);
    theta3=atan2(S3,C3)*180/pi;
    
    C2=((l1+l2*C3)*(pX(a)*C1+pY(a)*S1)-l2*pZ(a)*S3)/(pZ(a)^2+(pX(a)*C1+pY(a)*S1)^2);
    S2=(1-C2^2)^(1/2);
    theta2=atan2(S2,C2)*180/pi;
    
 
    T01=[cosd(theta1)  ,-sind(theta1) ,0  ,0;
         sind(theta1)  ,cosd(theta1)  ,0  ,0;
         0             ,0             ,1  ,0;
         0             ,0             ,0  ,1];
 
    T12=[cosd(theta2)  ,-sind(theta2) ,0  ,0;
         0             ,0             ,-1 ,0;
         sind(theta2)  ,cosd(theta2)  ,0  ,0;
         0             ,0             ,0  ,1];
 
    T23=[cosd(theta3)  ,-sind(theta3) ,0  ,l1;
         -sind(theta3) ,-cosd(theta3) ,0  ,0;
         0             ,0             ,-1 ,0;
         0             ,0             ,0  ,1];
 
    T34=[1   ,0   ,0  ,l2;
         0   ,1   ,0  ,0;
         0   ,0   ,1  ,0;
         0   ,0   ,0  ,1];
 
    T02=T01*T12;
    T03=T02*T23;
    T04=T03*T34;
 
    Lx1=[T02(1,4) T03(1,4)];
    Ly1=[T02(2,4) T03(2,4)];
    Lz1=[T02(3,4) T03(3,4)];
    
    Lx2=[T03(1,4) T04(1,4)];
    Ly2=[T03(2,4) T04(2,4)];
    Lz2=[T03(3,4) T04(3,4)];
    
    plot3(Lx1,Ly1,Lz1,'r','LineWidth',1)    %Dibujo del primer eslabòn
    grid on
    hold on
    axis([-20 40 -20 20 -0 20])
    plot3(Lx2,Ly2,Lz2,'b','LineWidth',1)    %Dibujo del segundo eslabòn
    hold on
    plot3(px,py,pz,'k','LineWidth',2)
    hold on
    plot3(pX,pY,pZ,'k','LineWidth',2)       %Dibujo de la recta eslabòn
    hold off
    
   
    pause(0.01)
end 
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%LINEA     L


%Programa de un robot manipuñador de GDL
%Trayectoria Circular parametrizada por polinomio de 5to grado
%Denavit-Hartenberg
%Dr. Waldemar Pérez Bailón

%%clc
%%clear all

l1=20;
l2=20;
 
p1=[0 15 0];
p2=[0 0 0];
 
x=p2(1)-p1(1);
y=p2(2)-p1(2);
z=p2(3)-p1(3);
 
largo=(x^2+y^2+z^2)^(1/2);
tf=50;
qf=largo;
t=0:0.1:tf;
q=qf*(10*(t.^3)/(tf^3)-15*(t.^4)/(tf^4)+6*(t.^5)/(tf^5));
 
 
for a=1:1:length(t)
    
    pXX(a)=p1(1)+x*q(a)/largo;
    pYY(a)=p1(2)+y*q(a)/largo;
    pZZ(a)=p1(3)+z*q(a)/largo;
    
    theta1=atan2(pYY(a),pXX(a))*180/pi;
    C1=cosd(theta1);
    S1=sind(theta1);
    
    C3=((pXX(a)^2)+(pYY(a)^2)+(pZZ(a)^2)-(l2^2)-(l1^2))/(2*l1*l2);
    S3=(1-C3^2)^(1/2);
    theta3=atan2(S3,C3)*180/pi;
    
    C2=((l1+l2*C3)*(pXX(a)*C1+pYY(a)*S1)-l2*pZZ(a)*S3)/(pZZ(a)^2+(pXX(a)*C1+pYY(a)*S1)^2);
    S2=(1-C2^2)^(1/2);
    theta2=atan2(S2,C2)*180/pi;
    
 
    T01=[cosd(theta1)  ,-sind(theta1) ,0  ,0;
         sind(theta1)  ,cosd(theta1)  ,0  ,0;
         0             ,0             ,1  ,0;
         0             ,0             ,0  ,1];
 
    T12=[cosd(theta2)  ,-sind(theta2) ,0  ,0;
         0             ,0             ,-1 ,0;
         sind(theta2)  ,cosd(theta2)  ,0  ,0;
         0             ,0             ,0  ,1];
 
    T23=[cosd(theta3)  ,-sind(theta3) ,0  ,l1;
         -sind(theta3) ,-cosd(theta3) ,0  ,0;
         0             ,0             ,-1 ,0;
         0             ,0             ,0  ,1];
 
    T34=[1   ,0   ,0  ,l2;
         0   ,1   ,0  ,0;
         0   ,0   ,1  ,0;
         0   ,0   ,0  ,1];
 
    T02=T01*T12;
    T03=T02*T23;
    T04=T03*T34;
 
    Lx1=[T02(1,4) T03(1,4)];
    Ly1=[T02(2,4) T03(2,4)];
    Lz1=[T02(3,4) T03(3,4)];
    
    Lx2=[T03(1,4) T04(1,4)];
    Ly2=[T03(2,4) T04(2,4)];
    Lz2=[T03(3,4) T04(3,4)];
    
    plot3(Lx1,Ly1,Lz1,'r','LineWidth',1)    %Dibujo del primer eslabòn
    grid on
    hold on
    axis([-20 40 -20 20 -0 20])
    plot3(Lx2,Ly2,Lz2,'b','LineWidth',1)    %Dibujo del segundo eslabòn
    hold on
    plot3(px,py,pz,'k','LineWidth',2)
    hold on
    plot3(pX,pY,pZ,'k','LineWidth',2)       %Dibujo de la recta eslabòn
    hold on
    plot3(pXX,pYY,pZZ,'k','LineWidth',2)       %Dibujo de la recta eslabòn
    hold off
    pause(0.01)
end 
%%%%%%%%SEGUNDO CIRCULO DE LA O
  %%%%%%%%%%%%%%%circulo de la a
l1=20;
l2=20;
cx=10;
cy=5;
cz=0;
r=5;
largo=2*pi*r;
tf=20;
qf=largo;
 
t=0:0.1:tf;
q=qf*(10*(t.^3)/(tf^3)-15*(t.^4)/(tf^4)+6*(t.^5)/(tf^5));
 
for a=1:1:length(t)
    pxx(a)=r*sind(360*q(a)/qf)+cx;
    pyy(a)=r*cosd(360*q(a)/qf)+cy;
    pzz(a)=0*r*cosd(360*q(a)/qf)+cz;
    
    theta1=atan2(pyy(a),pxx(a))*180/pi;
    C1=cosd(theta1);
    S1=sind(theta1);
    
    C3=((pxx(a)^2)+(pyy(a)^2)+(pzz(a)^2)-(l2^2)-(l1^2))/(2*l1*l2);
    S3=(1-C3^2)^(1/2);
    theta3=atan2(S3,C3)*180/pi;
    
    C2=((l1+l2*C3)*(pxx(a)*C1+pyy(a)*S1)-l2*pzz(a)*S3)/(pzz(a)^2+(pxx(a)*C1+pyy(a)*S1)^2);
    S2=(1-C2^2)^(1/2);
    theta2=atan2(S2,C2)*180/pi;
 
    T01=[cosd(theta1)  ,-sind(theta1) ,0  ,0;
         sind(theta1)  ,cosd(theta1)  ,0  ,0;
         0             ,0             ,1  ,0;
         0             ,0             ,0  ,1];
 
    T12=[cosd(theta2)  ,-sind(theta2) ,0  ,0;
         0             ,0             ,-1 ,0;
         sind(theta2)  ,cosd(theta2)  ,0  ,0;
         0             ,0             ,0  ,1];
 
    T23=[cosd(theta3)  ,-sind(theta3) ,0  ,l1;
         -sind(theta3) ,-cosd(theta3) ,0  ,0;
         0             ,0             ,-1 ,0;
         0             ,0             ,0  ,1];
 
    T34=[1   ,0   ,0  ,l2;
         0   ,1   ,0  ,0;
         0   ,0   ,1  ,0;
         0   ,0   ,0  ,1];
 
    T02=T01*T12;
    T03=T02*T23;
    T04=T03*T34;
 
    Lx1=[T02(1,4) T03(1,4)];
    Ly1=[T02(2,4) T03(2,4)];
    Lz1=[T02(3,4) T03(3,4)];
    Lx2=[T03(1,4) T04(1,4)];
    Ly2=[T03(2,4) T04(2,4)];
    Lz2=[T03(3,4) T04(3,4)];
    
    plot3(Lx1,Ly1,Lz1,'r','LineWidth',1)
    grid on
    hold on
    axis([-20 40 -20 20 -0 20])
    plot3(Lx2,Ly2,Lz2,'b','LineWidth',1)
    hold on
    plot3(px,py,pz,'k','LineWidth',2)
    hold on
    plot3(pX,pY,pZ,'k','LineWidth',2)       %Dibujo de la recta eslabòn
    hold on
    plot3(pXX,pYY,pZZ,'k','LineWidth',2)       %Dibujo de la recta eslabòn
    hold on
    plot3(pxx,pyy,pzz,'k','LineWidth',2)
    hold off
    pause(0.1)
end
%%%%%%%%%%%%%%%%%%palito de la segunda a
%%%%%%%%%%%%%%%%LINEA


%Programa de un robot manipuñador de GDL
%Trayectoria Circular parametrizada por polinomio de 5to grado
%Denavit-Hartenberg
%Dr. Waldemar Pérez Bailón

%%clc
%%clear all

l1=20;
l2=20;
 
p1=[15 10 0];
p2=[15 0 0];
 
x=p2(1)-p1(1);
y=p2(2)-p1(2);
z=p2(3)-p1(3);
 
largo=(x^2+y^2+z^2)^(1/2);
tf=50;
qf=largo;
t=0:0.1:tf;
q=qf*(10*(t.^3)/(tf^3)-15*(t.^4)/(tf^4)+6*(t.^5)/(tf^5));
 
 
for a=1:1:length(t)
    
    pXXX(a)=p1(1)+x*q(a)/largo;
    pYYY(a)=p1(2)+y*q(a)/largo;
    pZZZ(a)=p1(3)+z*q(a)/largo;
    
    theta1=atan2(pYYY(a),pXXX(a))*180/pi;
    C1=cosd(theta1);
    S1=sind(theta1);
    
    C3=((pXXX(a)^2)+(pYYY(a)^2)+(pZZZ(a)^2)-(l2^2)-(l1^2))/(2*l1*l2);
    S3=(1-C3^2)^(1/2);
    theta3=atan2(S3,C3)*180/pi;
    
    C2=((l1+l2*C3)*(pXXX(a)*C1+pYYY(a)*S1)-l2*pZZZ(a)*S3)/(pZZZ(a)^2+(pXXX(a)*C1+pYYY(a)*S1)^2);
    S2=(1-C2^2)^(1/2);
    theta2=atan2(S2,C2)*180/pi;
    
 
    T01=[cosd(theta1)  ,-sind(theta1) ,0  ,0;
         sind(theta1)  ,cosd(theta1)  ,0  ,0;
         0             ,0             ,1  ,0;
         0             ,0             ,0  ,1];
 
    T12=[cosd(theta2)  ,-sind(theta2) ,0  ,0;
         0             ,0             ,-1 ,0;
         sind(theta2)  ,cosd(theta2)  ,0  ,0;
         0             ,0             ,0  ,1];
 
    T23=[cosd(theta3)  ,-sind(theta3) ,0  ,l1;
         -sind(theta3) ,-cosd(theta3) ,0  ,0;
         0             ,0             ,-1 ,0;
         0             ,0             ,0  ,1];
 
    T34=[1   ,0   ,0  ,l2;
         0   ,1   ,0  ,0;
         0   ,0   ,1  ,0;
         0   ,0   ,0  ,1];
 
    T02=T01*T12;
    T03=T02*T23;
    T04=T03*T34;
 
    Lx1=[T02(1,4) T03(1,4)];
    Ly1=[T02(2,4) T03(2,4)];
    Lz1=[T02(3,4) T03(3,4)];
    
    Lx2=[T03(1,4) T04(1,4)];
    Ly2=[T03(2,4) T04(2,4)];
    Lz2=[T03(3,4) T04(3,4)];
    
    plot3(Lx1,Ly1,Lz1,'r','LineWidth',1)    %Dibujo del primer eslabòn
    grid on
    hold on
    axis([-20 40 -20 20 -0 20])
    plot3(Lx2,Ly2,Lz2,'b','LineWidth',1)    %Dibujo del segundo eslabòn
    hold on
    plot3(px,py,pz,'k','LineWidth',2)
    hold on
    plot3(pX,pY,pZ,'k','LineWidth',2)       %Dibujo de la recta eslabòn
    hold on
    plot3(pXX,pYY,pZZ,'k','LineWidth',2)       %Dibujo de la recta eslabòn
    hold on
    plot3(pxx,pyy,pzz,'k','LineWidth',2)
    hold on
    plot3(pXXX,pYYY,pZZZ,'k','LineWidth',2)
    hold off
    
    
    pause(0.01)
end


%%%%%%%%%%%%%%%%%%%%%%%%%PRIMER PALO DE LA N
l1=20;
l2=20;
 
p1=[20 10 0];
p2=[20 0 0];
 
x=p2(1)-p1(1);
y=p2(2)-p1(2);
z=p2(3)-p1(3);
 
largo=(x^2+y^2+z^2)^(1/2);
tf=50;
qf=largo;
t=0:0.1:tf;
q=qf*(10*(t.^3)/(tf^3)-15*(t.^4)/(tf^4)+6*(t.^5)/(tf^5));
 
 
for a=1:1:length(t)
    
    pXXXX(a)=p1(1)+x*q(a)/largo;
    pYYYY(a)=p1(2)+y*q(a)/largo;
    pZZZZ(a)=p1(3)+z*q(a)/largo;
    
    theta1=atan2(pYYYY(a),pXXXX(a))*180/pi;
    C1=cosd(theta1);
    S1=sind(theta1);
    
    C3=((pXXXX(a)^2)+(pYYYY(a)^2)+(pZZZZ(a)^2)-(l2^2)-(l1^2))/(2*l1*l2);
    S3=(1-C3^2)^(1/2);
    theta3=atan2(S3,C3)*180/pi;
    
    C2=((l1+l2*C3)*(pXXXX(a)*C1+pYYYY(a)*S1)-l2*pZZZZ(a)*S3)/(pZZZZ(a)^2+(pXXXX(a)*C1+pYYYY(a)*S1)^2);
    S2=(1-C2^2)^(1/2);
    theta2=atan2(S2,C2)*180/pi;
    
 
    T01=[cosd(theta1)  ,-sind(theta1) ,0  ,0;
         sind(theta1)  ,cosd(theta1)  ,0  ,0;
         0             ,0             ,1  ,0;
         0             ,0             ,0  ,1];
 
    T12=[cosd(theta2)  ,-sind(theta2) ,0  ,0;
         0             ,0             ,-1 ,0;
         sind(theta2)  ,cosd(theta2)  ,0  ,0;
         0             ,0             ,0  ,1];
 
    T23=[cosd(theta3)  ,-sind(theta3) ,0  ,l1;
         -sind(theta3) ,-cosd(theta3) ,0  ,0;
         0             ,0             ,-1 ,0;
         0             ,0             ,0  ,1];
 
    T34=[1   ,0   ,0  ,l2;
         0   ,1   ,0  ,0;
         0   ,0   ,1  ,0;
         0   ,0   ,0  ,1];
 
    T02=T01*T12;
    T03=T02*T23;
    T04=T03*T34;
 
    Lx1=[T02(1,4) T03(1,4)];
    Ly1=[T02(2,4) T03(2,4)];
    Lz1=[T02(3,4) T03(3,4)];
    
    Lx2=[T03(1,4) T04(1,4)];
    Ly2=[T03(2,4) T04(2,4)];
    Lz2=[T03(3,4) T04(3,4)];
    
    plot3(Lx1,Ly1,Lz1,'r','LineWidth',1)    %Dibujo del primer eslabòn
    grid on
    hold on
    axis([-20 40 -20 20 -0 20])
    plot3(Lx2,Ly2,Lz2,'b','LineWidth',1)    %Dibujo del segundo eslabòn
    hold on
    plot3(px,py,pz,'k','LineWidth',2)
    hold on
    plot3(pX,pY,pZ,'k','LineWidth',2)       %Dibujo de la recta eslabòn
    hold on
    plot3(pXX,pYY,pZZ,'k','LineWidth',2)       %Dibujo de la recta eslabòn
    hold on
    plot3(pxx,pyy,pzz,'k','LineWidth',2)
    hold on
    plot3(pXXX,pYYY,pZZZ,'k','LineWidth',2)
    hold on
    plot3(pXXXX,pYYYY,pZZZZ,'k','LineWidth',2)
    hold off
    
    
    pause(0.01)
end
  %%%%%%%%%%%%%%% circulo de la n
l1=20;
l2=20;
cx=25;
cy=5;
cz=0;
r=5;
largo=2*pi*r;
tf=20;
qf=largo;
 
t=0:0.1:tf;
q=qf*(10*(t.^3)/(tf^3)-15*(t.^4)/(tf^4)+6*(t.^5)/(tf^5));
 
for a=1:1:length(t)
    pxxx(a)=r*sind((180*q(a)/qf)-90)+cx;
    pyyy(a)=r*cosd((180*q(a)/qf)-90)+cy;
    pzzz(a)=0*r*cosd((180*q(a)/qf)-90)+cz;
    
    theta1=atan2(pyyy(a),pxxx(a))*180/pi;
    C1=cosd(theta1);
    S1=sind(theta1);
    
    C3=((pxxx(a)^2)+(pyyy(a)^2)+(pzzz(a)^2)-(l2^2)-(l1^2))/(2*l1*l2);
    S3=(1-C3^2)^(1/2);
    theta3=atan2(S3,C3)*180/pi;
    
    C2=((l1+l2*C3)*(pxxx(a)*C1+pyyy(a)*S1)-l2*pzzz(a)*S3)/(pzzz(a)^2+(pxxx(a)*C1+pyyy(a)*S1)^2);
    S2=(1-C2^2)^(1/2);
    theta2=atan2(S2,C2)*180/pi;
 
    T01=[cosd(theta1)  ,-sind(theta1) ,0  ,0;
         sind(theta1)  ,cosd(theta1)  ,0  ,0;
         0             ,0             ,1  ,0;
         0             ,0             ,0  ,1];
 
    T12=[cosd(theta2)  ,-sind(theta2) ,0  ,0;
         0             ,0             ,-1 ,0;
         sind(theta2)  ,cosd(theta2)  ,0  ,0;
         0             ,0             ,0  ,1];
 
    T23=[cosd(theta3)  ,-sind(theta3) ,0  ,l1;
         -sind(theta3) ,-cosd(theta3) ,0  ,0;
         0             ,0             ,-1 ,0;
         0             ,0             ,0  ,1];
 
    T34=[1   ,0   ,0  ,l2;
         0   ,1   ,0  ,0;
         0   ,0   ,1  ,0;
         0   ,0   ,0  ,1];
 
    T02=T01*T12;
    T03=T02*T23;
    T04=T03*T34;
 
    Lx1=[T02(1,4) T03(1,4)];
    Ly1=[T02(2,4) T03(2,4)];
    Lz1=[T02(3,4) T03(3,4)];
    Lx2=[T03(1,4) T04(1,4)];
    Ly2=[T03(2,4) T04(2,4)];
    Lz2=[T03(3,4) T04(3,4)];
    
    plot3(Lx1,Ly1,Lz1,'r','LineWidth',1)
    grid on
    hold on
    axis([-20 40 -20 20 -0 20])
    plot3(Lx2,Ly2,Lz2,'b','LineWidth',1)
    hold on
    plot3(px,py,pz,'k','LineWidth',2)
    hold on
    plot3(pX,pY,pZ,'k','LineWidth',2)       %Dibujo de la recta eslabòn
    hold on
    plot3(pXX,pYY,pZZ,'k','LineWidth',2)       %Dibujo de la recta eslabòn
    hold on
    plot3(pxx,pyy,pzz,'k','LineWidth',2)
    hold on
    plot3(pXXX,pYYY,pZZZ,'k','LineWidth',2)
    hold on
    plot3(pXXXX,pYYYY,pZZZZ,'k','LineWidth',2)
    hold on
    plot3(pxxx,pyyy,pzzz,'k','LineWidth',2)
    hold off
    pause(0.1)
end
%%%%%%%%%%%%%%%%%%%%%%%%%PRIMER PALO DE LA N
l1=20;
l2=20;
 
p1=[30 5 0];
p2=[30 0 0];
 
x=p2(1)-p1(1);
y=p2(2)-p1(2);
z=p2(3)-p1(3);
 
largo=(x^2+y^2+z^2)^(1/2);
tf=50;
qf=largo;
t=0:0.1:tf;
q=qf*(10*(t.^3)/(tf^3)-15*(t.^4)/(tf^4)+6*(t.^5)/(tf^5));
 
 
for a=1:1:length(t)
    
    pXXXXX(a)=p1(1)+x*q(a)/largo;
    pYYYYY(a)=p1(2)+y*q(a)/largo;
    pZZZZZ(a)=p1(3)+z*q(a)/largo;
    
    theta1=atan2(pYYYYY(a),pXXXXX(a))*180/pi;
    C1=cosd(theta1);
    S1=sind(theta1);
    
    C3=((pXXXXX(a)^2)+(pYYYYY(a)^2)+(pZZZZZ(a)^2)-(l2^2)-(l1^2))/(2*l1*l2);
    S3=(1-C3^2)^(1/2);
    theta3=atan2(S3,C3)*180/pi;
    
    C2=((l1+l2*C3)*(pXXXXX(a)*C1+pYYYYY(a)*S1)-l2*pZZZZZ(a)*S3)/(pZZZZZ(a)^2+(pXXXXX(a)*C1+pYYYYY(a)*S1)^2);
    S2=(1-C2^2)^(1/2);
    theta2=atan2(S2,C2)*180/pi;
    
 
    T01=[cosd(theta1)  ,-sind(theta1) ,0  ,0;
         sind(theta1)  ,cosd(theta1)  ,0  ,0;
         0             ,0             ,1  ,0;
         0             ,0             ,0  ,1];
 
    T12=[cosd(theta2)  ,-sind(theta2) ,0  ,0;
         0             ,0             ,-1 ,0;
         sind(theta2)  ,cosd(theta2)  ,0  ,0;
         0             ,0             ,0  ,1];
 
    T23=[cosd(theta3)  ,-sind(theta3) ,0  ,l1;
         -sind(theta3) ,-cosd(theta3) ,0  ,0;
         0             ,0             ,-1 ,0;
         0             ,0             ,0  ,1];
 
    T34=[1   ,0   ,0  ,l2;
         0   ,1   ,0  ,0;
         0   ,0   ,1  ,0;
         0   ,0   ,0  ,1];
 
    T02=T01*T12;
    T03=T02*T23;
    T04=T03*T34;
 
    Lx1=[T02(1,4) T03(1,4)];
    Ly1=[T02(2,4) T03(2,4)];
    Lz1=[T02(3,4) T03(3,4)];
    
    Lx2=[T03(1,4) T04(1,4)];
    Ly2=[T03(2,4) T04(2,4)];
    Lz2=[T03(3,4) T04(3,4)];
    
    plot3(Lx1,Ly1,Lz1,'r','LineWidth',1)    %Dibujo del primer eslabòn
    grid on
    hold on
    axis([-20 40 -20 20 -0 20])
    plot3(Lx2,Ly2,Lz2,'b','LineWidth',1)    %Dibujo del segundo eslabòn
    hold on
    plot3(px,py,pz,'k','LineWidth',2)
    hold on
    plot3(pX,pY,pZ,'k','LineWidth',2)       %Dibujo de la recta eslabòn
    hold on
    plot3(pXX,pYY,pZZ,'k','LineWidth',2)       %Dibujo de la recta eslabòn
    hold on
    plot3(pxx,pyy,pzz,'k','LineWidth',2)
    hold on
    plot3(pXXX,pYYY,pZZZ,'k','LineWidth',2)
    hold on
    plot3(pXXXX,pYYYY,pZZZZ,'k','LineWidth',2)
    hold on
    plot3(pxxx,pyyy,pzzz,'k','LineWidth',2)
    hold on
    plot3(pXXXXX,pYYYYY,pZZZZZ,'k','LineWidth',2)
    
    hold off
    
    
    pause(0.01)
end

