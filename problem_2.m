close all;
clear; clc;

%% heat equation: u_t = u_xx + u_yy
%x,y,t step; start and finish parameters
Jxy=50;
length_x=Jxy+1;
length_y=Jxy+1;
delta_xy = 1/Jxy;  
delta_t = 0.0008; 

x0=0;
xl=1;
y0=0;
yl=1;
t0=0;
tl=0.2;

x = x0:delta_xy:xl; 
y = y0:delta_xy:yl;
t = t0:delta_t:tl;

%make the whole plate into a super_plate in line
u=zeros(length_x*length_y,1);

% hot point
    %right & up select 4 point 
u(((length_y-8)*length_x+length_y-8):((length_y-8)*length_x+length_y-7))=1;
u(((length_y-9)*length_x+length_y-8):((length_y-9)*length_x+length_y-7))=1;
    %left & down select 4 point
u(7*length_x+8:7*length_x+9)=1;
u(8*length_x+8:8*length_x+9)=1;

% mark point in/out circle in specific matrix
% use the unequation: distance of(any _point to center) ≤ 1/3
mark_matrix = ones(length(x), length(y));
for j = 5:length(y)-5
    for i = 5:length(x)-5
        if ((j-1)*delta_xy-0.5)^2 +((i-1)*delta_xy-0.5)^2 <0.33^2
            mark_matrix(i,j) = 0; 
        end
    end
end

%% ADI & boundary process
u_time = Heat_Process(x,y,t,u,mark_matrix);

%% plot the picture 
t_arry_display = [1,7,26,63];
Draw_3D_Temperture(u_time,x,y,t_arry_display,delta_t,mark_matrix);

%% Heat Process
function u_time = Heat_Process(x,y,t,u,mark_matrix)

    % Initial the parameter miu,dt,dx,dy
    Jxy = 50;
    delta_xy = x(2) - x(1);
    delta_t = t(2) - t(1);
    miu = delta_t/(delta_xy^2);
    grid_mesh = Jxy+1; %grid length
    
    tl = 0.2;
    u_time(:,1) = u; 
    uc = zeros(Jxy+1,Jxy+1,length(t)+1);
    u_adi = zeros(Jxy+1,Jxy+1,length(t)+1);
    
    for n = 2:length(t)
        u_temp_pre = reshape(u_time(:,n-1),grid_mesh,grid_mesh); 
        u_temp_now = u_temp_pre; 
        for i = 2:grid_mesh-1 
            for j = 2:grid_mesh-1 
                if(mark_matrix(i,j)) 
                    if i+j>51 % x+y>1
                         %  case 1
                         if mark_matrix(i,j+1)==0 && mark_matrix(i-1,j)==0  
                            
                            lamda = Weight_Lamda_Calculate(i,j);
                            gamma = Weight_Gamma_Calculate(i,j);
                            
                            u_temp_now(i,j) = 2*miu*u_temp_pre(i+1,j)/(1+lamda) +2*miu*u_temp_pre(i,j-1)/(1+gamma)+...
                                (1-2*miu/gamma-2*miu/lamda)*u_temp_pre(i,j);
                            
                          % case 4
                          elseif mark_matrix(i,j-1)==0 && mark_matrix(i-1,j)==0 %4

                              lamda = Weight_Lamda_Calculate(i,j);
                              gamma = Weight_Gamma_Calculate(i,j);
                                                       
                              u_temp_now(i,j) = 2*miu/(1+gamma)*(u_temp_pre(i+1,j)-(1+1/gamma)*u_temp_pre(i,j))+ ...
                                  2*miu/(1+lamda)*(u_temp_pre(i,j+1)-(1+1/lamda)*u_temp_pre(i,j)) + u_temp_pre(i,j);
                          
                          % case 5  
                          elseif mark_matrix(i,j-1)==0 &&mark_matrix(i+1,j)==0 

                              lamda = Weight_Lamda_Calculate(i,j);
                              gamma = Weight_Gamma_Calculate(i,j);
                                                        
                              u_temp_now(i,j) = 2*miu/(1+gamma)*(u_temp_pre(i-1,j)-(1+1/gamma)*u_temp_pre(i,j))+ ...
                                2*miu/(1+lamda)*(u_temp_pre(i,j+1)-(1+1/lamda)*u_temp_pre(i,j)) + u_temp_pre(i,j);
  
                          %  case 6    
                          elseif mark_matrix(i,j-1)==0 &&mark_matrix(i+1,j) &&mark_matrix(i-1,j)
                            
                              lamda = Weight_Lamda_Calculate(i,j);
                                                        
                              u_temp_now(i,j) = miu*(u_temp_pre(i+1,j) - 2*u_temp_pre(i,j) + u_temp_pre(i-1,j))+ ...
                                2*miu/(1+lamda)*(u_temp_pre(i,j+1)-(1+1/lamda)*u_temp_pre(i,j)) + u_temp_pre(i,j);
                        
                          % case 7  
                          elseif mark_matrix(i-1,j)==0 && mark_matrix(i,j+1)&& mark_matrix(i,j-1)%7
                            
                            gamma = Weight_Gamma_Calculate(i,j);
                            
                            u_temp_now(i,j) = 2*miu/(1+gamma)*(u_temp_pre(i+1,j)-(1+1/gamma)*u_temp_pre(i,j))+ ...
                                miu*(u_temp_pre(i,j-1)-2*u_temp_pre(i,j)+u_temp_pre(i,j+1)) + u_temp_pre(i,j);
                            
                            
                         else %point not on the boundary of the surface and outside the circle
                             uc = ADI_y(u_adi,length(x)-1,delta_t,tl,u_temp_now,u_temp_pre,n,uc,j);
                         end
                         
%                       Neumann boundary condition
%                       as gradient is 0,  the former equation change as follow:
%                       gradient(B) =  (u_B-u_Z) /q*delta_y = 0  have u_B = u_Z
%                       gradient(D) =  (u_D-u_Z) /q*delta_x = 0  have u_D = u_Z
%                       horizon:    u_Z = p*u_W + (1-p)*u_P     p = |PZ|/delta_x
%                       vertical:   u_Z = p*u_N + (1-p)*u_P     p = |PZ|/delta_y	
                    else 
                        %   case 2
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

                        %   case 1
                        elseif mark_matrix(i,j+1)==0 &&mark_matrix(i-1,j)==0 
                            x_v_inter = (i-1)*delta_xy;                        
                            y_v_inter = 1/2 - sqrt(0.33*0.33 - (x_v_inter-1/2)^2);                       
                            lamda = abs(y_v_inter-(j-1)*delta_xy)/delta_xy;
                            y_left= (j-1)*delta_xy;                        
                            x_left = 1/2 + sqrt(0.33*0.33 - (y_left-1/2)^2);
                            gamma = abs(x_left-(i-1)*delta_xy)/delta_xy;
                          
                            z_v_x = (x_v_inter-1/2)/(y_v_inter-1/2)*((j-1)*delta_xy -1/2) + 1/2;
                            if (z_v_x-(i-1)*delta_xy) <= delta_xy % horizon
                         
                                p = abs(z_v_x-(i-1)*delta_xy)/delta_xy;
                                u_z_v = p*u_temp_pre(i+1,j)+(1-p)*u_temp_pre(i,j);
                            else % vertical 
                                z_v_y = (y_v_inter-1/2)/(x_v_inter-1/2)*((i)*delta_xy -1/2) + 1/2;
                                z_v_x = (i)*delta_xy;
                                p = abs(z_v_y - (j-1)*delta_xy)/delta_xy;
                                u_z_v = p*u_temp_pre(i+1,j+1)+(1-p)*u_temp_pre(i+1,j);
                            end
                            
                            z_h_x = (x_h_inter-1/2)/(y_h_inter-1/2)*((j-2)*delta_xy -1/2) + 1/2;
                            if (-z_h_x+ (i-1)*delta_xy) <= delta_xy % horizon
                                p = abs(z_v_x-(i-1)*delta_xy)/delta_xy;
                                u_z_h = p*u_temp_pre(i-1,j-1)+(1-p)*u_temp_pre(i,j-1);
                            else % vertical
                                z_h_y = (y_v_inter-1/2)/(x_v_inter-1/2)*((i-1)*delta_xy -1/2) + 1/2;
                                p = abs(z_h_y - (j-1)*delta_xy)/delta_xy;
                                u_z_h = p*u_temp_pre(i,j-1)+(1-p)*u_temp_pre(i,j);
                            end
                                                  
                            u_temp_now(i,j) = 2*miu/(1+gamma)*(u_temp_pre(i+1,j)-(1+1/gamma)*u_temp_pre(i,j)+1/gamma*u_z_h )+ ...
                                  2*miu/(1+lamda)*(u_temp_pre(i,j-1)-(1+1/lamda)*u_temp_pre(i,j)+1/lamda*u_z_v) + u_temp_pre(i,j);
                              
                        %   case 3      
                        elseif mark_matrix(i,j+1)==0 &&mark_matrix(i-1,j)&&mark_matrix(i+1,j)
                            x_v_inter = (i-1)*delta_xy;                        
                            y_v_inter = 1/2 - sqrt(0.33*0.33 - (x_v_inter-1/2)^2);                       
                            lamda = abs(y_v_inter-(j-1)*delta_xy)/delta_xy;
                          
                            z_v_x = (x_v_inter-1/2)/(y_v_inter-1/2)*((j-1)*delta_xy -1/2) + 1/2;
                            if (-z_v_x+(i-1)*delta_xy) <= delta_xy &&(-z_v_x+(i-1)*delta_xy)>0%The intersection of the grid is on the horizontal line
           
                                p = abs(z_v_x-(i-1)*delta_xy)/delta_xy;
                                u_z_v = p*u_temp_pre(i-1,j)+(1-p)*u_temp_pre(i,j);
                            elseif (-z_v_x+(i-1)*delta_xy)>delta_xy
                                z_v_y = (y_v_inter-1/2)/(x_v_inter-1/2)*((i-2)*delta_xy -1/2) + 1/2;
                               
                                p = abs(z_v_y - (j-1)*delta_xy)/delta_xy;
                                u_z_v = p*u_temp_pre(i-1,j+1)+(1-p)*u_temp_pre(i-1,j);
                            elseif (z_v_x-(i-1)*delta_xy) <= delta_xy && (z_v_x-(i-1)*delta_xy)>0
                                p = abs(z_v_x-(i-1)*delta_xy)/delta_xy;
                                u_z_v = p*u_temp_pre(i+1,j)+(1-p)*u_temp_pre(i,j);
                            elseif (z_v_x-(i-1)*delta_xy) > delta_xy %The intersection of the grid is on the vertical line
                                z_v_y = (y_v_inter-1/2)/(x_v_inter-1/2)*((i)*delta_xy -1/2) + 1/2;
                                
                                p = abs(z_v_y - (j-1)*delta_xy)/delta_xy;
                                u_z_v = p*u_temp_pre(i+1,j+1)+(1-p)*u_temp_pre(i+1,j);
                            end
                                  
                            u_temp_now(i,j) = miu*(u_temp_pre(i+1,j) - 2*u_temp_pre(i,j) + u_temp_pre(i-1,j))+ ...
                                      2*miu/(1+lamda)*(u_temp_pre(i,j-1)-(1+1/lamda)*u_temp_pre(i,j)+1/lamda*u_z_v) + u_temp_pre(i,j);
                                  
                        %   case 5
                        elseif mark_matrix(i,j-1)==0 && mark_matrix(i+1,j)==0 
                            x_down = (i-1)*delta_xy;                        
                            y_down = 1/2 + sqrt(0.33*0.33 - (x_down-1/2)^2);
                            lamda = abs(y_down-(j-1)*delta_xy)/delta_xy;
                            y_h_inter= (j-1)*delta_xy;                        
                            x_h_inter = 1/2 - sqrt(0.33*0.33 - (y_h_inter-1/2)^2);
                            gamma = abs(x_h_inter -(i-1)*delta_xy)/delta_xy;
                            
                            z_v_x = (x_down-1/2)/(y_down-1/2)*((j-1)*delta_xy -1/2) + 1/2;
                            if (-z_v_x+(i-1)*delta_xy) <= delta_xy % horizon
                                
                                p = abs(z_v_x-(i-1)*delta_xy)/delta_xy;
                                u_z_v = p*u_temp_pre(i-1,j)+(1-p)*u_temp_pre(i,j);
                            else % vertical 
                                z_v_y = (y_down-1/2)/(x_down-1/2)*((i-2)*delta_xy -1/2) + 1/2;
                               
                                p = abs(z_v_y - (j-1)*delta_xy)/delta_xy;
                                u_z_v = p*u_temp_pre(i-1,j-1)+(1-p)*u_temp_pre(i-1,j);
                            end
                            
                            z_h_y = (y_h_inter-1/2)/(x_h_inter-1/2)*((i-1)*delta_xy -1/2) + 1/2;
                            if (z_h_y- (j-1)*delta_xy) <= delta_xy % horizon
                                
                                p = abs(z_h_y- (j-1)*delta_xy)/delta_xy;
                                u_z_h = p*u_temp_pre(i,j+1)+(1-p)*u_temp_pre(i,j);
                            else % vertical 
                                z_h_x = (x_h_inter-1/2)/(y_h_inter-1/2)*((j)*delta_xy -1/2) + 1/2;
                                
                                p = abs(z_h_x - (i-1)*delta_xy)/delta_xy;
                                u_z_h = p*u_temp_pre(i+1,j+1)+(1-p)*u_temp_pre(i,j+1);
                            end
                            
                            u_temp_now(i,j) = 2*miu/(1+gamma)*(u_temp_pre(i-1,j)-(1+1/gamma)*u_temp_pre(i,j)+1/gamma*u_z_h )+ ...
                                  2*miu/(1+lamda)*(u_temp_pre(i,j+1)-(1+1/lamda)*u_temp_pre(i,j)+1/lamda*u_z_v) + u_temp_pre(i,j);

                        %   case 8
                        elseif mark_matrix(i+1,j)==0 &&mark_matrix(i,j+1)&&mark_matrix(i,j-1) 
                            y_h_inter= (j-1)*delta_xy;                        
                            x_h_inter = 1/2 - sqrt(0.33*0.33 - (y_h_inter-1/2)^2);
                            gamma = abs(x_h_inter -(i-1)*delta_xy)/delta_xy;
                            
                            z_h_y = (y_h_inter-1/2)/(x_h_inter-1/2)*((i-1)*delta_xy -1/2) + 1/2;
                            if (z_h_y- (j-1)*delta_xy) <= delta_xy &&(z_h_y- (j-1)*delta_xy) > 0 %horizon
                                
                                p = abs(z_h_y-(j-1)*delta_xy)/delta_xy;
                                u_z_h = p*u_temp_pre(i,j+1)+(1-p)*u_temp_pre(i,j);
                            elseif (z_h_y- (j-1)*delta_xy) > delta_xy
                                z_h_x = (x_h_inter-1/2)/(y_h_inter-1/2)*((j)*delta_xy -1/2) + 1/2;
                                
                                p = abs(z_h_x-(i-1)*delta_xy)/delta_xy;
                                u_z_h =  p*u_temp_pre(i+1,j+1)+(1-p)*u_temp_pre(i,j+1);
                            elseif -z_h_y+ (j-1)*delta_xy<= delta_xy && -z_h_y+ (j-1)*delta_xy>0
                                
                                p = abs(z_h_y-(j-1)*delta_xy)/delta_xy;
                                u_z_h = p*u_temp_pre(i,j-1)+(1-p)*u_temp_pre(i,j);
                            elseif -z_h_y+ (j-1)*delta_xy> delta_xy
                                z_h_x = (x_h_inter-1/2)/(y_h_inter-1/2)*((j-2)*delta_xy -1/2) + 1/2;
                                
                                p = abs(z_h_x-(i-1)*delta_xy)/delta_xy;
                                u_z_h =  p*u_temp_pre(i+1,j-1)+(1-p)*u_temp_pre(i,j-1);
                            end
                                                    
                            u_temp_now(i,j) = 2*miu/(1+gamma)*(u_temp_pre(i-1,j)-(1+1/gamma)*u_temp_pre(i,j)+1/gamma*u_z_h)+ ...
                                miu*(u_temp_pre(i,j-1)-2*u_temp_pre(i,j)+u_temp_pre(i,j+1)) + u_temp_pre(i,j);
                                                   
                        else %point not on the boundary of the surface and outside the circle
                            uc = ADI_y(u_adi,length(x)-1,delta_t,tl,u_temp_now,u_temp_pre,n,uc,j);
                        end
                    end                    
                end
            end
            u_temp = ADI_x(u_adi,length(x)-1, delta_t, tl, u_temp_now, n, uc,i); 
            u_temp_now(:,:)= u_temp(:,:,n);
        end       
        % Neumann boundary
        %  line boundary 
        u_temp_now(:,1) = u_temp_now(:,2);
        u_temp_now(1,:) = u_temp_now(2,:);
        %  curve boundary 
        u_time(:,n) = u_temp_now(:);
    end
end

%% ADI-method
function uc = ADI_y(u_adi,Jxy,delta_t,tl,u_temp,u_temp_pre,k,uc,j)
%*******************************
%   calculate numerical solution in PR scheme
%   u_adi the numerical solution  
    %list related para
    delta_xy = 1/Jxy;
    miu = delta_t/delta_xy^2;
      
    % initial
    u_adi(:,:,k) = u_temp(:,:);
    u_adi(:,:,k-1) = u_temp_pre(:,:);
    A_calculate = diag(ones(Jxy-1,1)*(1+miu),0)+diag(ones(Jxy-2,1)*(-miu/2),-1)+diag( ...
        ones(Jxy-2,1)*(-miu/2),1);

    % in y order，make k layer to k+1/2 layer
    % initial 
    % uc in pre
    a1 = miu*u_adi(1,j-1,k)-2*miu*u_adi(1,j,k)+miu*u_adi(1,j+1,k);
    b1 = miu*u_adi(1,j-1,k-1)-2*miu*u_adi(1,j,k-1)+miu*u_adi(1,j+1,k-1);
    uc(1,j,k) = 1/2*(u_adi(1,j,k-1)+u_adi(1,j,k))-1/4*(a1-b1);
    % uc in last
    am = miu*u_adi(Jxy+1,j-1,k)-2*miu*u_adi(Jxy,j-1,k)+miu*u_adi(Jxy,j+1,k);
    bm = miu*u_adi(Jxy+1,j-1,k-1)-2*miu*u_adi(Jxy,j-1,k-1)+miu*u_adi(Jxy,j+1,k-1);
    uc(Jxy,j,k)=1/2*(u_adi(Jxy,j,k-1)+u_adi(Jxy,j,k))-1/4*(am-bm);
    % list pre layer matrix to calculate middle layer --y
    B1 = miu/2*u_adi(2:Jxy,j-1,k-1)+(1-miu)*u_adi(2:Jxy,j,k-1)+miu/2*u_adi(2:Jxy,j+1,k-1);

    %calculate process, B1 has already uplevel to a specific matrix (notice dilemention of B1)
    B1(1) = B1(1)+miu/2*uc(1,j,k);
    B1(Jxy-1) = B1(Jxy-1) + miu/2*uc(Jxy+1,j,k);
   % B1 = B1';
    uc(2:Jxy,j,k) = A_calculate\B1;  %   1 & Jxy+1 has insert before calculate
end

function u_adi = ADI_x(u_adi,Jxy,delta_t,tl,u_temp,k,uc,i)
%*******************************
%   calculate numerical solution in PR scheme
%   u_adi the numerical solution 

    %list related para
    delta_xy = 1/Jxy;
    miu = delta_t/delta_xy^2;
    
    % initial
    u_adi(:,:,k) = u_temp(:,:);
    A_calculate = diag(ones(Jxy-1,1)*(1+miu),0)+diag(ones(Jxy-2,1)*(-miu/2),-1)+diag( ...
        ones(Jxy-2,1)*(-miu/2),1);

    % in x order，make k+1/2 layer to k+1 layer
    % list k+1/2 layer matrix to calculate bext layer --x
    B2 = miu/2*uc(i-1,2:Jxy,k)+(1-miu)*uc(i,2:Jxy,k)+miu/2*uc(i+1,2:Jxy,k);

    %calculate process, B2 has already uplevel to a specific matrix (notice dilemention of B2)
    B2(1) = B2(1)+miu/2*u_adi(i,1,k);
    B2(Jxy-1) = B2(Jxy-1)+miu/2*u_adi(i,Jxy+1,k);
    B2 = reshape(B2,[Jxy-1,1]); % equal to transpose
    u_adi(i,2:Jxy,k) = A_calculate\B2;
end

%% calculate weight
function lamda = Weight_Lamda_Calculate(i,j)
    delta_xy = 1/50;
    x_v_inter = (i-1)*delta_xy;
    y_v_inter = 1/2 - sqrt(0.33*0.33 - (x_v_inter-1/2)^2);
    lamda = abs(y_v_inter-(j-1)*delta_xy)/delta_xy;
end

function gamma = Weight_Gamma_Calculate(i,j)
    delta_xy = 1/50;
    y_h_inter = (j-1)*delta_xy;
    x_h_inter = 1/2 + sqrt(0.33*0.33 - (y_h_inter -1/2)^2);
    gamma = abs(x_h_inter-(i-1)*delta_xy)/delta_xy;
end

%% Draw temperture
function Draw_3D_Temperture(u_time, x,y,t_arry_display,delta_t,mark_matrix)

    for i = 1:length(t_arry_display) 
        mesh_temp_draw = reshape(u_time(:,t_arry_display(i)),length(x),length(y)); 
        for j = 5:length(y)-5
            for i_x = 5:length(x)-5
                if mark_matrix(i_x,j) == 0
                    mesh_temp_draw(i_x,j) = -0.1*max(mesh_temp_draw(:));
                end
            end
        end    
        % Output the maximum temperature at that moment
        disp(max(mesh_temp_draw(:)));
        t_arr = [0, 0.005, 0.02, 0.05];
        % plot figure
        subplot(2,2,i)
        s=mesh(x,y,mesh_temp_draw);
        colormap hot; 
        colorbar;
        s.FaceColor = 'flat';
        title("t = "+ t_arr(i));
        xlabel("x");
        ylabel("y");
        zlabel("Tem")
    end      
end