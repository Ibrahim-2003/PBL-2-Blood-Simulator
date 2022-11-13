% N = 200;
% T1 =  2;
% T2 = 4;
% T3 = 6;
% T4 = 7;
% alpha_max = 30;
% am = alpha_max;
% t = 0 :T4/(N-1):T4;
% % Implement the piecewise function
% y1 = (10/T1^3)*t.^3-(15/T1^4)*t.^4+(6/T1^5)*t.^5;
% y1 = am*y1;
% 
% y2 = am;
% 
% y3 = (10/T1^3)*(T3-t).^3-(15/T1^4)*(T3-t).^4+(6/T1^5)*(T3-t).^5;
% y3 = am*y3;
% 
% y4 = 0;
% alpha = y1.*(0<t & t<T1) +y2.*(T1<t & t<T2)+y3.*(T2<t & t<T3)+...
%         y4.*(T1<t & t<T2);
% % Implement the derivative
% y5 = (30/T1^3)*t.^2-(60/T1^4)*t.^3+(30/T1^5)*t.^4;
% y5 = am*y5;
% y6 = 0;
% y7 = (-30/T1^3)*(T3-t).^2+(60/T1^4)*(T3-t).^3-(30/T1^5)*(T3-t).^4;
% y7 = am*y7;
% y8 = 0;
% alpha_prime = y5.*(0<t & t<T1) +y6.*(T1<t & t<T2)+y7.*(T2<t & t<T3)+...
%         y8.*(T1<t & t<T2);  
% figure
% plot(t ,alpha,'linewidth',2)
% hold on
% plot(t ,alpha_prime,'linewidth',2,'color','r')
% legend('\alpha','\alpha_{prime}');
% a = title('\alpha(t)');
% set(a,'fontsize',14);
% a = ylabel('y');
% set(a,'Fontsize',14);
% a = xlabel('t [0 7]');
% set(a,'Fontsize',14);
% xlim([0 T4])
% grid

Twosecond_oder_odes();

function Twosecond_oder_odes 
% Verify that your equations are implemented correctly.
% Run these functions as a single m-file. 
N = 200;
T = 12;
t = 0 :T/(N-1):T;
ic = [.5 .5 1 0]; % Initial conditions Must be 4 numbers, try different values
[t,z]= ode45(@RHS, t,ic );

figure
plot(t ,z(:,1),'linewidth',2,'color','b')
hold on
plot(t ,z(:,2),'linewidth',2,'color','r')
plot(t ,z(:,3),'linewidth',2,'color','m')
plot(t ,z(:,4),'linewidth',2,'color','g')
legend('\beta','\beta_{prime}','x',',x_{prime}');
a = title('\beta(t) x(t)');
set(a,'fontsize',14);
a = ylabel('y');
set(a,'Fontsize',14);
a = xlabel('t [0 7]');
set(a,'Fontsize',14);
ylim([-18000 16000])
grid

end

function dZdt= RHS(t,z)
% parameters
J2 = 0.90;
J3 = 0.10;
mA = 10;
mB = 35;
r1 = 0.25;
r2 = 0.55;
r3 = 0.15;
cc = 300;
c  =  400;
kc = 2000;
k  = 6000;
T1 =  2;
T2 = 4;
T3 = 6;
alpha_max = 30;
am = alpha_max;

% Coefficients from Equation 1
a1 = J2 + 2*J3 + mA*r3^2;
a2 = 2*cc*r2^2 + c*r3^2;
a3 = -c*r3;
a4 = 2*kc*r2^2 + k*r3^2;
a5 = -k*r3;
% Coefficients from Equation 2
b1 = mB;
b2 = -c*r3;
b3 = c;
b4 = a5;
b5 = k;
% Implements the piecewise function, alpha (Test this to make sure its
% right) 
y1 = (10/T1^3)*t.^3-(15/T1^4)*t.^4+(6/T1^5)*t.^5;
y1 = am*y1;
y2 = am;
y3 = (10/T1^3)*(T3-t).^3-(15/T1^4)*(T3-t).^4+(6/T1^5)*(T3-t).^5;
y3 = am*y3;
y4 = 0;
alpha = y1.*(0<t & t<T1) +y2.*(T1<t & t<T2)+y3.*(T2<t & t<T3)+...
        y4.*(T1<t & t<T2);
% Implements the derivative of alpha (Do it manually!)
y5 = (30/T1^3)*t.^2-(60/T1^4)*t.^3+(30/T1^5)*t.^4;
y5 = am*y5;
y6 = 0;
y7 = (-30/T1^3)*(T3-t).^2+(60/T1^4)*(T3-t).^3-(30/T1^5)*(T3-t).^4;
y7 = am*y7;
y8 = 0;
alpha_prime = y5.*(0<t & t<T1) +y6.*(T1<t & t<T2)+y7.*(T2<t & t<T3)+...
        y8.*(T1<t & t<T2);

   % This is the right hand side of Equation 1 
     F = 2*cc*r1*r2*alpha_prime + 2*kc*r1*r2*alpha;

     % Implements two 2nd orders into 4 1st orders ODES
     dZdt_1 = z(1); % Solution for beta (z1 = beta, z2 = derivative of beta)
     dZdt_3 = z(3); % Solution for x (z3 = x, z4 = derivative of x)

     dZdt_2 = a2*z(2)+a4*z(1)+a3*z(4)+a5*z(3)-F; % dZdt_2 = double derivative of beta
     dZdt_2 = (-1/a1)*dZdt_2;

     dZdt_4 = b2*z(2)+b4*z(1)+b3*z(4)+b5*z(3);% dZdt_4 = double derivative of x
     dZdt_4 = (-1/b1)*dZdt_4;
     dZdt =[dZdt_1; dZdt_2; dZdt_3;dZdt_4];

end