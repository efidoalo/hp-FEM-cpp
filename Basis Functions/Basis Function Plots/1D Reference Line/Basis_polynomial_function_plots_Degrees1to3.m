x=0:0.05:1;
B1Linear=zeros(21);
B2Linear=zeros(21);
B1Quadratic=zeros(21);
B2Quadratic=zeros(21);
B3Quadratic=zeros(21);
B1Cubic=zeros(21);
B2Cubic=zeros(21);
B3Cubic=zeros(21);
B4Cubic=zeros(21);
for i=1:21
    B1Linear(i) = 1 - x(i);
    B2Linear(i) = x(i);
    B1Quadratic(i) = 2*x(i)*x(i) - 3*x(i) + 1;
    B2Quadratic(i) = -4*x(i)*x(i) + 4*x(i);
    B3Quadratic(i) = 2*x(i)*x(i) - x(i);
    B1Cubic(i) = -4.5*x(i)*x(i)*x(i) + 9*x(i)*x(i) -5.5*x(i) + 1;
    B2Cubic(i) = 13.5*x(i)*x(i)*x(i) -22.5*x(i)*x(i) + 9*x(i);
    B3Cubic(i) = -13.5*x(i)*x(i)*x(i) + 18*x(i)*x(i) - 4.5*x(i);
    B4Cubic(i) = 4.5*x(i)*x(i)*x(i) - 4.5*x(i)*x(i) + x(i);
end

figure(1);
plot(x,B1Linear,x,B2Linear);
title('Basis function polynomials of the finite element FE_{1}=([0,1],1)');
xlabel('x');
ylabel('y');

figure(2);
plot(x,B1Quadratic, x,B2Quadratic, x, B3Quadratic);
title('Basis function polynomials of the finite element FE_{2}=([0,1],2)');
xlabel('x');
ylabel('y');

figure(3);
plot(x,B1Cubic, x,B2Cubic, x, B3Cubic, x, B4Cubic);
title('Basis function polynomials of the finite element FE_{3}=([0,1],3)');
xlabel('x');
ylabel('y');

figure(4)
subplot(2,1,1);
plot(x,B1Linear,x,B2Linear);
text11='\phi_{k2}'; text12='\phi_{k1}';
text(0.2,0.3,text11); text(0.76,0.3,text12);
title('Basis function polynomials of the finite element FE_{1}=([0,1],1)');
xlabel('x');
subplot(2,1,2);
plot(x,B1Quadratic, x,B2Quadratic, x, B3Quadratic);
text21='\phi_{k1}'; text22='\phi_{k2}'; text23='\phi_{k3}';
text(0.075,0.83,text21); text(0.48,0.9,text22); text(0.88,0.83,text23);
title('Basis function polynomials of the finite element FE_{2}=([0,1],2)');
xlabel('x');
axis([0 1 -0.2 1.1]);


