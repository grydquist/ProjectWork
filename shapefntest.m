x=[1 -3 2;5 2 1;1 0 5;-1 -2 -3]';
N=[1 0 0;0 1 0;0 0 1;-1 -1 -1]';
xi=zeros(3,3);
for a=1:4
    xi(:,1) = xi(:,1)+x(:,a)*N(1,a);
    xi(:,2) = xi(:,2)+x(:,a)*N(2,a);
    xi(:,3) = xi(:,3)+x(:,a)*N(3,a);
end

xinv=inv(xi);

xco=[.2;-1;-1];
diff=xco-x(:,4);
shps=xinv*diff

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
