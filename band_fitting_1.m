clear;
eps1=0.605;
eps2=1.972;
t0=-0.169;
t1=0.228;
t2=0.390;
t11=0.207;
t12=0.239;
t22=0.252;

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