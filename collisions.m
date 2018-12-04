d1 = .4;
d2 = 1.6;
k=1;
m1 = .6;
m2 = .6;
x1=[1,2,0];
x2=[1,1,0];
v1=[1,-2,0];
v2=[1,3,0];

n1=(x1 - x2)/(d1+d2)*2;
n2=-n1;


t1=cross(cross(n1,v1),n1)/norm(cross(cross(n1,v1),n1));
t2=cross(cross(n2,v2),n2)/norm(cross(cross(n2,v2),n2));

vperp1=dot(t1,v1);
vperp2=dot(t2,v2);
vpar1 =dot(v1,n1);
vpar2 =dot(v2,n2);

pa=m1*vpar1-m2*vpar2;
pb=(-vpar1 - vpar2)*k;

upar2 = (pa-m1*pb)/(m1+m2);
upar1=pb+upar2;
upar2=-upar2;

vf1 = vperp1*t1 + upar1*n1
vf2 = vperp2*t2 + upar2*n2

plot([-v1(1),0],[-v1(2),0])
hold on
plot([0,n1(1)],[0,n1(2)])
plot([vf1(1),0],[vf1(2),0])


