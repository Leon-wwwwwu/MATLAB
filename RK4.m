%���ݽṹ��1000������/����100s
x = zeros(1,1001);  %�߶�
y = zeros(1,1001);  %�ٶ�
z = zeros(1,1001);  %����
t = zeros(1,1001); %ʱ��
m = zeros(1,1001);
%��ʼֵ
i=2;
m(1)=9;
h    = 0.1;         %����
t(1) = 0;           %��ʼʱ��
x(1) = 80000;       %��ʼ��ظ߶�
y(1) = -5000;       %��ʼ�ٶ�
z(1) = 200000;      %��ʼ����
g0   = 9.8;         %��ʼ�������ٶ�
r0   = 6371000;     %����뾶
T0   = 7600000;     %�����������������
Isp  = 342;         %�ȳ�
f1   = @(x, y, z, t)(y); %����1
f2   = @(x, y, z, t)(-1*g0*(r0/(r0+x))^2 +num(t)*T0/z + fk(y)*y*abs(y)/z);%����2
f3   = @(x, y, z, t)(-1*num(t)*T0/(Isp*g0));%����3
%��������
for  i = 2 : 1001
    t(i) = roundn((t(i-1) + h),-1);%����ʱ��
    m(i)=num(t(i));
    k11 =  f1(x(i-1),y(i-1),z(i-1),t(i-1));
    k21 =  f2(x(i-1),y(i-1),z(i-1),t(i-1));
    k31 =  f3(x(i-1),y(i-1),z(i-1),t(i-1));
    
    k12 =  f1(x(i-1)+0.5*h*k11,y(i-1)+0.5*h*k21,z(i-1)+0.5*h*k31,t(i-1)+0.5*h);
    k22 =  f2(x(i-1)+0.5*h*k11,y(i-1)+0.5*h*k21,z(i-1)+0.5*h*k31,t(i-1)+0.5*h);
    k32 =  f3(x(i-1)+0.5*h*k11,y(i-1)+0.5*h*k21,z(i-1)+0.5*h*k31,t(i-1)+0.5*h);
    
    k13 =  f1(x(i-1)+0.5*h*k12,y(i-1)+0.5*h*k22,z(i-1)+0.5*h*k32,t(i-1)+0.5*h);
    k23 =  f2(x(i-1)+0.5*h*k12,y(i-1)+0.5*h*k22,z(i-1)+0.5*h*k32,t(i-1)+0.5*h);
    k33 =  f3(x(i-1)+0.5*h*k12,y(i-1)+0.5*h*k22,z(i-1)+0.5*h*k32,t(i-1)+0.5*h);
    
    k14 =  f1(x(i-1)+h*k13,y(i-1)+h*k23,z(i-1)+h*k33,t(i-1)+h);
    k24 =  f2(x(i-1)+h*k13,y(i-1)+h*k23,z(i-1)+h*k33,t(i-1)+h);
    k34 =  f3(x(i-1)+h*k13,y(i-1)+h*k23,z(i-1)+h*k33,t(i-1)+h);

    x(i) = x(i-1)+h*(1/6*k11 +1/3*k12+1/3*k13+1/6*k14);%���¸߶�
    y(i) = y(i-1)+h*(1/6*k21 +1/3*k22+1/3*k23+1/6*k24);%�����ٶ�
    z(i) = z(i-1)+h*(1/6*k31 +1/3*k32+1/3*k33+1/6*k34);%��������        
end
 figure(1);
 plot(t,x);
 xlabel('t'), ylabel('r','FontName','Times New Roman','FontSize',14,'Rotation',0);
 title('f(t, r)����');
 figure(2);
 plot(t,y);
 xlabel('t'), ylabel('v/(m/s)','FontName','Times New Roman','FontSize',14,'Rotation',0);
 title('f(t, v)����');
 figure(3);
 plot(t,z);
 xlabel('t'), ylabel('m/(kg)','FontName','Times New Roman','FontSize',14,'Rotation',0);
 title('f(t, m)����');

function n = num(t1)
    if t1 < 7.4
        n = 9;
    elseif ((t1>=68.0) &&(t1<69.4))||((t1>=73.4) && (t1<74.0))
        n = 2;
    elseif (t1>=74.6) && (t1<74.9)
        n=1; 
    else
        n = 0;
    end
end


function [k1] = fk(y)
  vx=[0 500 2000 5000];
  k=[1e-5 1.5e-5 1.8e-5 2e-5];
  if y<-5000
      k1=2e-5;
  else
      k1=interp1(vx,k,-y);
  end
end




