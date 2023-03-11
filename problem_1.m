clear
clc
%u(0,t)=u(1,t)=0
%u(x,0)=2x if 0<=x<=1/2 2-2x if1/2<=x<=1
%ut = uxx
Jx=100;
Jy=100;
delta_x=1/Jx;
delta_y=1/Jy;
delta_t=0.00002;

r=1;
s=1;

t0=0;
tl=0.1;

x0=0;
xl=1;
x=x0;
y0=0;
yl=1;
y=y0;

u=zeros(xl/delta_x+1,yl/delta_y+1,1000+1);
miu_x=delta_t/(delta_x.^2);
miu_y=delta_t/(delta_y.^2);

for n=1:tl/delta_t
    t=(n-1)*delta_t;     
   for r= 1:(xl/delta_x)
            x=(r-1)*delta_x;
        if r==1||r==Jx
                u(r,s,n)=0;        
        else
            for s=1:(yl/delta_y)
               y=(s-1)*delta_y;
                if n==1
                    if ((x^2+y^2)>=0.6&&(x^2+y^2)<=0.8)
                        u(r,s,n)=1;
                    end 
                else
                    if s==1||s==Jy
                        u(r,s,n)=0;
                    else
                        u(r,s,n)=miu_x*(u(r+1,s,n-1)-2*u(r,s,n-1)+u(r-1,s,n-1))+miu_x*(u(r,s+1,n-1)-2*u(r,s,n-1)+u(r,s-1,n-1))+u(r,s,n-1);    
                    end
                end  
            end
        end
    end
end
data_x=0:delta_x:1;
data_y=0:delta_y:1;
subplot(2,2,2);
surf(data_x,data_y,u(:,:,51));
title("﻿t = 0.001 for the data ");
xlabel("x");
ylabel("y");
zlabel("Temperture");
colorbar;
disp(max(u(:,:,51)));

subplot(2,2,3);
surf(data_x,data_y,u(:,:,201));
colorbar;
title("﻿t = 0.004 for the data ");
xlabel("x");
ylabel("y");
zlabel("Temperture");
disp(max(u(:,:,201)));

subplot(2,2,4);
surf(data_x,data_y,u(:,:,501));
colorbar;
title("﻿t = 0.01 for the data ");
xlabel("x");
ylabel("y");
zlabel("Temperture");
disp(max(u(:,:,501)));

subplot(2,2,1)
surf(data_x,data_y,u(:,:,1));
colorbar;
title("t = 0 ﻿ Initial data");
xlabel("x");
ylabel("y");
zlabel("Temperture");
disp(max(u(:,:,1)));
