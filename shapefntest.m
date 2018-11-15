x=[-0.5,-1.01,6.33;-0.51 -1.53 5.82;-1.93 -0.49 6.24;-0.55 -1.49 6.37]';
N=[1 0 0;0 1 0;0 0 1;-1 -1 -1]';
xi=zeros(3,3);
dummyx=zeros(4,4);
dummyb=dummyx;

Jac=zeros(4,4);
Jac(:,1)=1;
Jac(:,2:4)=x';
jd=det(Jac);
dummyx(:,1)=x(1,:)';
dummyx(:,2:4)=x';
a=1/jd*det(dummyx);


for a=1:4
    xi(:,1) = xi(:,1)+x(:,a)*N(1,a);
    xi(:,2) = xi(:,2)+x(:,a)*N(2,a);
    xi(:,3) = xi(:,3)+x(:,a)*N(3,a);
end

xinv=inv(xi);

xco=[-.6;-1.4;6.35];
diff=xco-x(:,4);
prnts=xinv*diff;
shps=zeros(1,4);
shps(1:3)=prnts;
shps(4)=1-sum(prnts);
shps

xp=zeros(1,4);
yp=xp;
zp=yp;
for a=1:4
    xp(a)=x(1,a);
    yp(a)=x(2,a);
    zp(a)=x(3,a);
end

T=[1,2,3;1,2,4;2,3,4;1,3,4];
trimesh(T,xp,yp,zp)
hold on
scatter3(xco(1),xco(2),xco(3))
