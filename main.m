for L=0.5: 0.5: 5
alpha=-0.6-0.3*1i;
tol=0.1;
h=0.01;
err=1;
while err>tol
f=rk(alpha,L);
f1=rk(alpha+h,L);
f2=rk(alpha+h*1i,L);
%define jacobian
J11=real((f1-f)/h);
J12=real((f2-f)/h);
J21=imag((f1-f)/h);
J22=imag((f2-f)/h);
J=[J11 J12; J21 J22];
A=inv(J);
F=[real(f); imag(f)];
alpha1=[real(alpha); imag(alpha)];
alpha2=alpha1-A*F;
disp(alpha2);
alpha=alpha2(1,1) + alpha2(1,1)*1i;
err=(alpha2(1,1)-alpha1(1,1))*(alpha2(1,1)-alpha1(1,1)) + (alpha2(2,1)-alpha1(2,1))*(alpha2(2,1)-alpha1(2,1));
err=sqrt(err);
alpha1=alpha2;
disp(err);
end
plot(L,imag(alpha),'xr');
hold on
end
