clear;
eps1=0.1;
eps2=0.1;
t0=0.1;
t1=0.401;
t2=0.1;
t11=0.1;
t12=0.1;
t22=0.1;

deps1=0.01;
deps2=0.01;
dt0=0.01;
dt1=0.01;
dt2=0.01;
dt11=0.01;
dt12=0.01;
dt22=0.01;

E1_gam_i=0;
E2_gam_i=2.9;
E3_gam_i=2.9;
E1_K_i=0;
E2_K_i=1.6;
E3_K_i=3.5;
E1_M_i=0;
E2_M_i=2.1;
E3_M_i=2.6;

for i=1:100000
    E_gam_1=eps1+6*t0;
    E_gam_2=eps2+3*(t11+t22);
    E_gam_3=eps2+3*(t11+t22);

    E_K_1=eps2-1.5*(t11+t22)-3^1.5*t12;
    E_K_2=eps1-3*t0;
    E_K_3=eps2-1.5*(t11+t22)+3^1.5*t12;

    f1=0.5*(eps1+eps2)-t0-1.5*t11+0.5*t22;
    f2=0.5*sqrt((eps1-eps2-2*t0+3*t11-t22)^2+64*t2^2);
    E_M_1=f1-f2;
    E_M_2=eps2+t11-3*t22;
    E_M_3=f1+f2;
    
    d=(E_gam_1-E1_gam_i)^2+(E_gam_2-E2_gam_i)^2+(E_gam_3-E3_gam_i)^2+(E_K_1-E1_K_i)^2+(E_K_2-E2_K_i)^2+(E_K_3-E3_K_i)^2+(E_M_1-E1_M_i)^2+(E_M_2-E2_M_i)^2+(E_M_3-E3_M_i)^2;
    dd_deps1=2*deps1*(E_gam_1-E1_gam_i)+2*deps1*(E_K_2-E2_K_i)+2*deps1*(E_M_1-E1_M_i)*(0.5-0.5*(eps1-eps2-2*t0+3*t11-t22)/(2*f2))+2*deps1*(E_M_3-E3_M_i)*(0.5+0.5*(eps1-eps2-2*t0+3*t11-t22)/(2*f2));
    dd_deps2=2*deps2*(E_gam_2-E2_gam_i)+2*deps2*(E_gam_3-E3_gam_i)+2*deps2*(E_K_1-E1_K_i)+2*deps2*(E_K_3-E3_K_i)+2*deps2*(E_M_1-E1_M_i)*(0.5-0.5*(eps1-eps2-2*t0+3*t11-t22)/(2*f2))+2*deps2*(E_M_2-E2_M_i)+2*deps2*(E_M_3-E3_M_i)*(0.5-0.5*(eps1-eps2-2*t0+3*t11-t22)/(2*f2));
    dd_dt0=2*dt0*(E_gam_1-E1_gam_i)+2*dt0*(E_K_3-E3_K_i)*(-3)+2*dt0*(E_M_1-E1_M_i)*(-1-0.5*(eps1-eps2-2*t0+3*t11-t22)*(-2)/(2*f2))+2*dt0*(E_M_3-E3_M_i)*(-1+0.5*(eps1-eps2-2*t0+3*t11-t22)*(-2)/(2*f2));
    dd_dt1=2*dt1*0;
    dd_dt2=-2*dt2*(E_M_1-E1_M_i)*(0.5*64*t2/(2*f2))+2*dt2*(E_M_3-E3_M_i)*(0.5*64*t2/(2*f2));
    dd_dt11=2*dt11*(E_gam_2-E2_gam_i)*3+2*dt11*(E_gam_3-E3_gam_i)*3+2*dt11*(E_K_1-E1_K_i)*(-1.5)+2*dt11*(E_K_3-E3_K_i)*(-1.5)+2*dt11*(E_M_1-E1_M_i)*(-1.5-0.5*(eps1-eps2-2*t0+3*t11-t22)*(3)/(2*f2))+2*dt11*(E_M_2-E2_M_i)+2*dt11*(E_M_3-E3_M_i)*(-1.5+0.5*(eps1-eps2-2*t0+3*t11-t22)*(3)/(2*f2));
    dd_dt12=-2*dt12*(E_K_1-E1_K_i)*(3^1.5)+2*dt12*(E_K_3-E3_K_i)*(3^1.5);
    dd_dt22=2*dt22*(E_gam_2-E2_gam_i)*3+2*dt22*(E_gam_3-E3_gam_i)*3+2*dt22*(E_K_1-E1_K_i)*(-1.5)+2*dt22*(E_K_3-E3_K_i)*(-1.5)+2*dt22*(E_M_1-E1_M_i)*(0.5-0.5*(eps1-eps2-2*t0+3*t11-t22)*(-1)/(2*f2))+2*dt22*(E_M_2-E2_M_i)*(-3)+2*dt22*(E_M_3-E3_M_i)*(0.5+0.5*(eps1-eps2-2*t0+3*t11-t22)*(-1)/(2*f2));
    
    eps1=eps1-dd_deps1;
    eps2=eps2-dd_deps2;
    t0=t0-dd_dt0;
    t1=t1-dd_dt1;
    t2=t2-dd_dt2;
    t11=t11-dd_dt11;
    t12=t12-dd_dt12;
    t22=t22-dd_dt22;
    disp(i)
end
disp(eps1)
disp(eps2)
disp(t0)
disp(t1)
disp(t2)
disp(t11)
disp(t12)
disp(t22)


for j=1:100
    alpha(j)=2*pi*j/(3*100);
    beta(j)=0;
    alpha(100+j)=2*pi/3-pi/6*(j/100);
    beta(j+100)=pi/2*(j/100);
    alpha(200+j)=pi/2-pi/2*(j/100);
    beta(200+j)=pi/2-pi/2*(j/100);
end

for k=1:300
    H(1,1)=2*t0*(cos(2*alpha(k))+2*cos(alpha(k))*cos(beta(k)))+eps1;
    H(1,2)=-2*sqrt(3)*t2*sin(alpha(k))*sin(beta(k))+2*1i*t1*(sin(2*alpha(k))+sin(alpha(k))*cos(beta(k)));
    H(1,3)=2*t2*(cos(2*alpha(k))-cos(alpha(k))*cos(beta(k)))+2*sqrt(3)*1i*t1*cos(alpha(k))*sin(beta(k));
    H(2,2)=2*t11*cos(2*alpha(k))+(t11+3*t22)*cos(alpha(k))*cos(beta(k))+eps2;
    H(3,3)=2*t22*cos(2*alpha(k))+(3*t11+t22)*cos(alpha(k))*cos(beta(k))+eps2;
    H(2,3)=sqrt(3)*(t22-t11)*sin(alpha(k))*sin(beta(k))+4*1i*t12*sin(alpha(k))*(cos(alpha(k))-cos(beta(k)));
    H(2,1)=conj(H(1,2));
    H(3,1)=conj(H(1,3));
    H(3,2)=conj(H(2,3));
    E=sort(eig(H));
    E1(k)=E(1);
    E2(k)=E(2);
    E3(k)=E(3);
end
xax=(1:1:300);
plot(xax(2:299),E1(2:299))
hold on
plot(xax(2:299),E2(2:299))
plot(xax(2:299),E3(2:299))