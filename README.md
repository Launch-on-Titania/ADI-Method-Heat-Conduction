# ADI_method_compute_heat_conduction
Both Neumann &amp; Dirichlet boundary condition are provided. Shown in 2 part and visualized by 3D mesh
# Problem 1
## Achieve Homogeneous Dirichlet Boundary in Unit Square

the eqution are as followed:

u(0,t)=u(1,t)=0

u(x,0)=2x if 0<=x<=1/2 2-2x if1/2<=x<=1

ut = uxx

At first, based on the the listed requirment.

make ``Jx=100;
Jy=100;
delta_x=1/Jx;
delta_y=1/Jy;
delta_t=0.00002;`` as initial parameters

then initial related matrix like``u=zeros(xl/delta_x+1,yl/delta_y+1,1000+1);``

and then, accomplish 3 times circulation(cause of explicit, there is no time limited):

1. discuss boundary;
```
if r==1||r==Jx
  u(r,s,n)=0;  
else % for circulation
  if n==1
    if ((x^2+y^2)>=0.6&&(x^2+y^2)<=0.8)
        u(r,s,n)=1;
    end 
  else
    if s==1||s==Jy
        u(r,s,n)=0;
```
2. the core explicit achieve:
   ```
     u(r,s,n)=miu_x*(u(r+1,s,n-1)-2*u(r,s,n-1)+u(r-1,s,n-1))+miu_x*(u(r,s+1,n-1)-2*u(r,s,n-1)+u(r,s-1,n-1))+u(r,s,n-1); 
   ```
finally, plot the surf picture in 4 situation.

<img src="/result/question1.bmp" width = "80%" /> 

# Problem 2 
## ADI Method for Heat Calculation

### processdure

At first, based on the the listed requirment.

make ``Jxy=50; length_x=Jxy+1; 
 length_y=Jxy+1; 
delta_xy = 1/Jxy;  
delta_t = 0.0008; `` as initial parameters

lists stripe of **x,y,t**,  ``x = x0:delta_xy:xl; 
y = y0:delta_xy:yl;
t = t0:delta_t:tl;``

then, initial beginning hot point
```
%right & up select 4 point 
u(((length_y-8)*length_x+length_y-8):((length_y-8)*length_x+length_y-7))=1;
u(((length_y-9)*length_x+length_y-8):((length_y-9)*length_x+length_y-7))=1;
%left & down select 4 point
u(7*length_x+8:7*length_x+9)=1;
u(8*length_x+8:8*length_x+9)=1;
```
 based on the center circle as null, marks point in/out circle in specific matrix:

 (use the unequation: distance of(any _point to center) ≤ 1/3)
 ```
mark_matrix = ones(length(x), length(y));
for j = 5:length(y)-5
    for i = 5:length(x)-5
        if ((j-1)*delta_xy-0.5)^2 +((i-1)*delta_xy-0.5)^2 <0.33^2
            mark_matrix(i,j) = 0; 
        end
    end
end
 ```
then, solve the temperature of the plate at each moment by CK ADI scheme method in normal domain:

here, the ADI method is nested in ``function u_time = Heat_Process(x,y,t,u,mark_matrix)``

after calculation of it, transports the related data to:

 ``function Draw_3D_Temperture(u_time, x,y, t_arry_display,delta_t,mark_matrix)``

 which achieve the final heat flow picture, creatively achieve colorbar in hot, which is more specific.

1. introduce the boundary of Dirichlet & Neumann boundary condition:
    
    1. the boundary for the circle divided in 8 case, the Dirichlet & Neumann owns 5 cases dependently.
        Dirichlet at first.
        
        the 2 conditions are distinguished by``i+j>51 % x+y>1 ``, by using mark_matrix sufficently, 
       achieves classfied. For example, in case 1:

            ```
            if mark_matrix(i,j+1)==0 && mark_matrix(i-1,j)==0                      
                lamda = Weight_Lamda_Calculate(i,j);
                gamma = Weight_Gamma_Calculate(i,j);
                u_temp_now(i,j) = 2*miu*u_temp_pre(i+1,j)/(1+lamda) +2*miu*u_temp_pre(i,j-1)/(1+gamma)+1-2*miu/gamma-2*miu/lamda)*u_temp_pre(i,j);
            end
            ```
    where the boundary is used explicit to process(which is limited quantity, for improving computing volecity)
    
    <img src="/result/heat_adi.png" width = "80%" /> 
      
    1. Neumann is talked as following:
          
       1. formula in progamm:
      
             Neumann boundary condition
                       as gradient is 0,  the former equation change as follow:

                       gradient(B) =  (u_B-u_Z) /q*delta_y = 0  have u_B = u_Z

                       gradient(D) =  (u_D-u_Z) /q*delta_x = 0  have u_D = u_Z

                       horizon:    u_Z = p*u_W + (1-p)*u_P     p = |PZ|/delta_x

                       vertical:   u_Z = p*u_N + (1-p)*u_P     p = |PZ|/delta_y	

         by using mark_matrix sufficently, 
       achieves classfied. For example, in case 2:

        ```
            if mark_matrix(i,j+1)==0 &&mark_matrix(i+1,j)==0 
                            x_v_inter = (i-1)*delta_xy;                        
                            y_v_inter = 1/2 - sqrt(0.33*0.33 - (x_v_inter-1/2)^2);                       
                            lamda = abs(y_v_inter-(j-1)*delta_xy)/delta_xy;
                            y_h_inter= (j-1)*delta_xy;                        
                            x_h_inter = 1/2 - sqrt(0.33*0.33 - (y_h_inter-1/2)^2);
                            gamma = abs(x_h_inter -(i-1)*delta_xy)/delta_xy;

                           z_v_x = (x_v_inter-1/2)/(y_v_inter-1/2)*((j-1)*delta_xy -1/2) + 1/2;
                            if (-z_v_x+(i-1)*delta_xy) <= delta_xy %  horizon

                                p = abs(z_v_x-(i-1)*delta_xy)/delta_xy;
                                u_z_v = p*u_temp_pre(i-1,j)+(1-p)*u_temp_pre(i,j);
                            else % vertical 
                                z_v_y = (y_v_inter-1/2)/(x_v_inter-1/2)*((i-2)*delta_xy -1/2) + 1/2;
                                z_v_x = (i-2)*delta_xy;
                                p = abs(z_v_y - (j-1)*delta_xy)/delta_xy;
                                u_z_v = p*u_temp_pre(i-1,j+1)+(1-p)*u_temp_pre(i-1,j);
                            end
                            
                            z_h_x = (x_h_inter-1/2)/(y_h_inter-1/2)*((j-2)*delta_xy -1/2) + 1/2;                           
                            if (z_h_x- (i-1)*delta_xy) <= delta_xy &&(z_h_x- (i-1)*delta_xy) > delta_xy % horizon
                                
                                p = abs(z_v_x-(i-1)*delta_xy)/delta_xy;
                                u_z_h = p*u_temp_pre(i+1,j-1)+(1-p)*u_temp_pre(i,j-1);
                            else % vertical 
                                z_h_y = (y_v_inter-1/2)/(x_v_inter-1/2)*((i-1)*delta_xy -1/2) + 1/2;
                                p = abs(z_h_y - (j-1)*delta_xy)/delta_xy;
                                u_z_h = p*u_temp_pre(i,j-1)+(1-p)*u_temp_pre(i,j);
                            end
                            
                            u_temp_now(i,j) = 2*miu/(1+gamma)*(u_temp_pre(i-1,j)-(1+1/gamma)*u_temp_pre(i,j)+1/gamma*u_z_h )+ ...
                                  2*miu/(1+lamda)*(u_temp_pre(i,j-1)-(1+1/lamda)*u_temp_pre(i,j)+1/lamda*u_z_v) + u_temp_pre(i,j);
       ```
2. on line boundary, use gradient = 0, 
   
   ``%  line boundary 
        u_temp_now(:,1) = u_temp_now(:,2);
        u_temp_now(1,:) = u_temp_now(2,:);`` 
3. once the mark matrix is beyond boundary, adi is called; ADI is in 2 steps;
   1. pre layer matrix to calculate middle layer --y``function uc = ADI_y(u_adi,Jxy,delta_t,tl,u_temp,u_temp_pre,k,uc,j)``
   2. k+1/2 layer matrix to calculate bext layer --x``function u_adi = ADI_x(u_adi,Jxy,delta_t,tl,u_temp,k,uc,i)``
   
   which is achieved one in y circulation, the other is in x circulation(calculate numerical solution in PR scheme).
