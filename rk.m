function d=rk(alpha,L)
%alpha=-0.1-0.1i;
Re=100;
k=1.0;
%calculate values of u1 and u2 matrix by rk4 method
a2=-(Re*alpha*alpha)/L -2*k*k;
a3=Re*alpha*alpha*k*k/L + k*k*k*k;
F1=@(t,y)[y(2);y(3);y(4);(-a2*y(3)-a3*y(1))];
[t,x]=ode45(F1,[-1 0],[0 0 1 0]);
u1=x(end,1:4);
%disp(u1);
F2=@(z,y)[y(2);y(3);y(4);(-a2*y(3)-a3*y(1))];
[z,sol]=ode45(F2,[-1 0],[0 0 0 1]);
u2=sol(end,1:4);
%disp(u2);

%calcultae values of v1 and v2 by rk4
%b2=-(Re/L)*((alpha+L*k*t2)*1i) -2*k*k;
%b3=(Re*k*k/L)*((alpha+L*k*t2)*1i) + k*k*k*k; 
F3=@(t2,y)[y(2);y(3);y(4);(-(-(Re/L)*((alpha+L*k*t2)*1i) -2*k*k)*y(3)-((Re*k*k/L)*((alpha+L*k*t2)*1i) + k*k*k*k)*y(1))];
[t2,x2]=ode45(F3,[1 0],[0 0 1 0]);
v1=x2(end,1:4);
%disp(v1);
F4=@(z2,y)[y(2);y(3);y(4);(-(-(Re/L)*((alpha+L*k*z2)*1i) -2*k*k)*y(3)-((Re*k*k/L)*((alpha+L*k*z2)*1i) + k*k*k*k)*y(1))];
[z2,sol2]=ode45(F4,[1 0],[0 0 0 1]);
v2=sol2(end,1:4);
%disp(v2);

% define ux, vx, pressures, sigma, tau to form detrerminant 
p1=-v1(1,2)*(alpha+k*k)/(k*k) + v1(1,4)+1i*(Re/k)*(v1(1,1));
p2=-v2(1,2)*(alpha+k*k)/(k*k) + v2(1,4)+1i*(Re/k)*(v2(1,1));

szz1= -p1+2*u1(1,2);
szz2= -p2+2*u2(1,2);
sxz1= (1i/k)*u1(1,3)+ 1i*k*u1(1,1);
sxz2= (1i/k)*u2(1,3)+ 1i*k*u2(1,1);
tauzz1 = -p1 + 2*v1(1,2);
tauzz2 = -p2 + 2*v2(1,2);
tauxz1 = 1i*k*v1(1,1) + (1i/k)*v1(1,3);
tauxz2 = 1i*k*v2(1,1) + (1i/k)*v2(1,3);
A21 = (1i*alpha/k)*u1(1,2)-L*u1(1,1);
A22 = (1i*alpha/k)*u2(1,2)-L*u2(1,1);
vx1 = (1i/k)*v1(1,2);
vx2 = (1i/k)*v2(1,2);
%dterminant
A=[alpha*u1(1,1) alpha*u2(1,1) v1(1,1) v2(1,1) ;A21 A22 vx1 vx2;sxz1 sxz2 tauxz1 tauxz2;szz1 szz2 tauzz1 tauzz2];
d=det(A);
%disp(d);
end
