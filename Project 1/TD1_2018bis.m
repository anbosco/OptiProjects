%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          MECA0027 Structural and Multidisciplinary Optimization         %
%                         Tutorial Class - 1                              %
%                     Unconstrained Optimization                          %
%                    University Of Liège, Belgium                         %
%                       Academic year 2018-2019                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Solve the problem
%             min f(x,y)= 2x^2 + 2y^2 - 3xy + 2y^2 - 2x + 10y - 1
%                                  or
%             min f(x,y)= 2x^4 + 2y^2 - 3xy + 2y^2 - 2x + 10y - 1
%
% using the following optimization methods

% 1) Steepest descent method
% 2) Conjugate direction method:
%   2.1) Conjugate gradient method
%   2.2) Using Fletcher-Reeves approximation
% 3) Newton method
% 4) Quasi-Newton method - Boyden-Fletcher-Goldfard-Shanno (BFGS)


function [OptSol,FVal] = TD1_2018bis(xinit)
close all
disp('Which optimization method do you want to use? Press:')
disp('     1 for Steepest descent method')
disp('     21 for Conjugate gradient method')
disp('     22 for Conjugate direction method using Fletcher-Reeves approximation')
disp('     3 for Newton method')
disp('     4 for Quasi-Newton method - Boyden-Fletcher-Goldfard-Shanno (BFGS)')
prompt = '';
ind = input(prompt);

% Parameters --------------------------------------------------------------------------
functionID = 2;             % =1 if min f(x,y)= 2x^2 + 2y^2 - 3xy + 2y^2 - 2x + 10y - 1
                            % =2 if min f(x,y)= 2x^4 + 2y^2 - 3xy + 2y^2 - 2x + 10y - 1
MaxIter = 100;              %Maximum number of iterations
n=2;                        %Dimension of the problem
xinit=reshape(xinit,2,1);   %Be sure it is a column vector
Epsilon=1e-4;               %Tolerance
x=zeros(n,MaxIter);         %Initialization of x
x(:,1)=xinit;               %Put xinit in vector x.
OF=zeros(1,MaxIter);         % Values of the objective function
nstep_bissection = zeros(1,MaxIter); % Number of step required by the bissection algorithm

% Methods ------------------------------------------------

if ind==1
    %% Steepest descent
    disp(' You chose the steepest descent method.');    
    for i=1:MaxIter
        if(i==1)  % We comppute the Hessian and the gradient at the current iterate
        H=getHess(x(:,i), functionID);
        df = getSens(x(:,i),functionID); 
        OF(i) = getObjFVal(x(:,i),functionID);
        end
        s = -df;        % We use the opposite of the gradient at the current iterate as descent direction
        [alpha,n_step] = getalpha(x(:,i),s,df,H,functionID);
        x(:,i+1) = x(:,i) +(alpha*s);       % Update of the iterate, the objective function value and its gradient  
        nstep_bissection(i+1) = n_step;
        OF(i+1) = getObjFVal(x(:,i+1),functionID);
        df = getSens(x(:,i+1),functionID);
        if(norm(df)<Epsilon)                % Verify if the next iterate is a stationnary point of f.  
            txt = sprintf('\n\n %s %s \n\n','Number of itérations : ', num2str(i));
            disp(txt);
            break;
        elseif(i==MaxIter)      
            txt = sprintf('\n\n %s %s \n\n','Impossible to reach the required precision in ', num2str(i));
            error(txt);
        end
        H=getHess(x(:,i+1), functionID);    % Update the Hessian
    end    
    x=x(:,1:i+1); %Remove the zero elements due to the initialization step    
    OF=OF(1,1:i+1);
    nstep_bissection = nstep_bissection(1,1:i+1);
    
     % Plot interesting information about the convergence
    if(functionID==1)
        for j = 1: i
           rho(j) = abs(OF(j+1)-OF(end))/(abs(OF(j)-OF(end)));
        end
         n_lin = log((abs(OF(j+1)+22.142857142857142)/(abs(OF(1)+22.142857142857142))))/log(rho(round(i/2)))
          n_iter_th = log((abs(OF(j+1)+22.142857142857142)/(abs(OF(1)+22.142857142857142))))/log(0.5625) % Only ok for f1
    else
        for j = 1: i
         rho(j) = abs(OF(j+1)-OF(end))/abs(OF(j)-OF(end));
        end
        n_lin = log((abs(OF(j)-OF(end))/(abs(OF(1)-OF(end)))))/log(rho(round(i/2)))
        
        Step_of_bissection = sum(nstep_bissection)
    end
    
elseif ind==21
    %% Gradient method (Fletcher and Reeves)
    disp(' You chose the conjugate gradient method.')
    for i=1:MaxIter
        if(i==1)    % Initialisation of the algorithm
            H = getHess(x(:,i), functionID);% Useful only for f1 (line search)
            df_k = getSens(x(:,i),functionID);
            d_k = -df_k;
        else        % Update of the descent direction (Fletcher and Reeves)
            num = norm(df_k)^2;
            denom = norm(df_k_1)^2;
            beta = num/denom;
            d_k = -df_k + beta*d_k_1;            
        end
        [alpha,n_step] = getalpha(x(:,i),d_k,df_k,H,functionID);    % Line search method
        x(:,i+1) = x(:,i) +(alpha*d_k);                             % Update of the iterate
        d_k_1 = d_k;                                                
        df_k_1 = df_k;
        df_k = getSens(x(:,i+1),functionID);                        % Update of the gradient     
        if(norm(df_k)<Epsilon)  
            txt = sprintf('\n\n %s %s \n\n','Number of itérations : ', num2str(i));
            disp(txt);
            break;
        elseif(i==MaxIter)      
            txt = sprintf('\n\n %s %s \n\n','Impossible to reach the required precision in ', num2str(i));
            error([txt]);
        end        
    end
    x=x(:,1:i+1); %Remove the zero elements due to the initialization step

elseif ind==22
    %% Gradient method (Polak-Ribiere)
    disp(' You chose the conjugate direction method using Fletcher-Reeves approximation.')
    
    for i=1:MaxIter
        if i==1 % Initialisation of the algorithm
            H=getHess(x(:,i), functionID);% Useful only for f1 (line search)
            df_k = getSens(x(:,i),functionID);   
            d_k = -df_k;
            OF(i) = getObjFVal(x(:,i),functionID);
        else % Update of the descent direction (Polak-Ribiere)
            num = dot(df_k,(df_k-df_k_1));
            denom = norm(df_k_1)^2;
            beta = num/denom;
            d_k = -df_k + beta*d_k_1; 
        end
        [alpha,n_step] = getalpha(x(:,i),d_k,df_k,H,functionID);    % Line search method
        x(:,i+1) = x(:,i) +(alpha*d_k);                             % Update of the iterate
        nstep_bissection(i+1) = n_step;
        OF(i+1) = getObjFVal(x(:,i+1),functionID);
        d_k_1 = d_k;
        df_k_1 = df_k;
        df_k = getSens(x(:,i+1),functionID);                        % Update of the gradient   
        if(norm(df_k)<Epsilon)
            txt = sprintf('\n\n %s %s \n\n','Number of itérations : ', num2str(i));
            disp(txt);
            break;
        elseif(i==MaxIter)      
            txt = sprintf('\n\n %s %s \n\n','Impossible to reach the required precision in ', num2str(i));
            error([txt]);
        end
    end
    x=x(:,1:i+1); %Remove the zero elements due to the initialization step
    OF=OF(1,1:i+1);
    nstep_bissection = nstep_bissection(1,1:i+1);

     % Plot interesting information about the convergence
    if(functionID==2)
        for j = 1: i
         err(j) = (abs(OF(j)-OF(end)));
         err_1(j) = abs(OF(j+1)-OF(end));
        end
        order = polyfit(log(err(1:end-1)),log(err_1(1:end-1)),1)      
        Step_of_bissection = sum(nstep_bissection)
    end
elseif ind==3
    %% Newton 
    disp(' You chose the Newton method.');
    for i=1:MaxIter
        mu = 0.001;
        if(i==1)
            H=getHess(x(:,i), functionID);
            df = getSens(x(:,i),functionID);    
        end
        vp = eig(H);
        while(vp(1)<=0||vp(2)<=0)
            H = H+mu*eye(2,2);
            vp = eig(H);
        end
        A = H^(-1);     
        alpha = 1;
        if(dot(df,-A*df)<0)  % Verify if it is a descent direction
            s = -A*df;
        else
            eig(H)
            s = A*df; 
        end
        while (getObjFVal(x(:,i) + alpha*s ,functionID) > getObjFVal(x(:,i),functionID))
            alpha = 0.5*alpha;
        end
        x(:,i+1) = x(:,i) + alpha*s;  
        H=getHess(x(:,i+1), functionID);
        df = getSens(x(:,i+1),functionID);
        if(norm(df)<Epsilon)
            txt = sprintf('\n\n %s %s \n\n','Number of itérations : ', num2str(i));
            disp(txt);
            break;
        elseif(i==MaxIter)      
            txt = sprintf('\n\n %s %s \n\n','Impossible to reach the required precision in ', num2str(i));
            error([txt]);
        end 
    end
    x=x(:,1:i+1); %Remove the zero elements due to the initialization step
    	
elseif ind==4
    %% BFGS
    disp(' You chose the quasi-Newton method (Boyden-Fletcher-Goldfard-Shanno (BFGS)).')
    for i=1:MaxIter
        if(i==1)
            H=eye(n);
            df = getSens(x(:,i),functionID);
        end  
        vp = eig(H);
        if(vp(1)<=0||vp(2)<=0)
           error('La theorie est fausse') 
        end
       s = -H*df;
       if functionID == 2
           [alpha,n_step] = getalpha(x(:,i),s,df,H,functionID);
       else
           [alpha,n_step] = getalpha(x(:,i),s,df,H,12);
       end
       x(:,i+1) = x(:,i) + alpha*s;
       delta = x(:,i+1)-x(:,i);       
       df_futur = getSens(x(:,i+1),functionID);
       gamma = df_futur-df;
       df = df_futur;
       num1 = gamma.'*H*gamma;
       denom1 = delta.'*gamma;
       num2 = delta*delta.';
       denom2 = delta.'*gamma;
       num3 = delta*gamma.'*H+H*gamma*delta.';
       denom3 = delta.'*gamma;
       H = H+(1+(num1/denom1))*(num2/denom2)-(num3/denom3);
       if(norm(df)<Epsilon)
            txt = sprintf('\n\n %s %s \n\n','Number of itérations : ', num2str(i));
            disp(txt);
            break;
%         elseif(i==MaxIter)      
%             txt = sprintf('\n\n %s %s \n\n','Impossible to reach the required precision in ', num2str(i));
%             error([txt]);
        end 
    end
    x=x(:,1:i+1); %Remove the zero elements due to the initialization step
     
else
    disp('Error: Your choice does not match the proposed methods.')
    return
end

OptSol=x(:,end);
disp(['The optimal point is: x = ' num2str(OptSol(1)) ', y = ' num2str(OptSol(2)) '.'])
FVal=getObjFVal(OptSol,functionID);
disp(['The objective function value is: ' num2str(FVal) '.'])

%% Plot the function with the optimization path.
if(functionID==2)
    lbx=-5;
    upx=4;
    lby=-7;
    upy=2;
else
    lbx=-12.5;
    upx=7.5;
    lby =lbx;
    upy = upx;
end
xi=lby:0.1:upy;
f=zeros(length(xi),length(xi));
for i=1:length(xi)
    for j=1:length(xi)
        f(j,i)=getObjFVal([xi(i);xi(j)],functionID);
    end
end
figure('name','Optimization path')
Figure1=figure(1);clf;set(Figure1,'defaulttextinterpreter','latex');
hold on;
set(gca,'fontsize',30,'fontname','Times','LineWidth',0.5);
hold on
axis equal
xlabel('$x_1$')
ylabel('$x_2$')
axis([lbx upx lby upy])
[C,h]=contour(xi,xi,f,[-10:2:10 10:10:50 75:25:300],'linewidth', 2);
colorbar;
%  clabel(C,h);
hold on
ind=1;
for i=1:size(x,2)-1
    plot(x(1,i),x(2,i),'.b','markersize',50)
    plot([x(1,i) x(1,i+1)],[x(2,i) x(2,i+1)],'-.k','linewidth',2)
    text(x(1,i),x(2,i),num2str(ind-1),'horizontalalignment','center','verticalalignment','middle','FontSize',14)
    ind=ind+1;
end
plot(x(1,end),x(2,end),'.r','markersize',50)
text(x(1,end),x(2,end),num2str(ind-1),'horizontalalignment','center','verticalalignment','middle','FontSize',14)

figure('name','Objective function')
Figure2=figure(2);clf;set(Figure2,'defaulttextinterpreter','latex');
hold on;
set(gca,'fontsize',30,'fontname','Times','LineWidth',0.5);
plot(0:size(transpose(OF))-1, transpose(OF),'r','linewidth',3);
plot(0:size(transpose(OF))-1,transpose(OF),'k.','markersize',25);
grid;
ylabel('Objective function');
xlabel('Number of iteration');

function fval = getObjFVal(x,functionID)
 if functionID == 1
     fval = 2*(x(1))^2-3*(x(1)*x(2))+2*(x(2))^2 - 2*x(1)+10*x(2)-1;
 elseif functionID == 2
     fval = 2*x(1)^4 - 3*x(1)*x(2) + 2*x(2)^2 - 2*x(1) + 10*x(2) - 1;
 end


 function df = getSens(x, functionID)
 if functionID == 1
     df(1) = 4*x(1) - 3*x(2) - 2;
     df(2) = -3*x(1) + 4*x(2) + 10;
     df = df.';
 elseif functionID == 2
    df(1) = 8*(x(1)^3)-3*x(2)-2;
    df(2) = -3*x(1)+4*x(2)+10;
    df = df.';
 end

 function H = getHess(x,functionID)
 if functionID == 1
     H(1,1) = 4;
     H(1,2) = -3;
     H(2,1) = -3;
     H(2,2) = 4;
 elseif functionID == 2
    H(1,1) = 24*(x(1)^2);
    H(1,2) = -3;
    H(2,1) = -3;
    H(2,2) = 4;
 end

 function [alpha,n_step] = getalpha(x_init,s,df,H,functionID)
     n_step = 0;
     %%%%%%%%%%%%%%%%%%
 if functionID == 1
     % The function 1 being a stricly convex quadratic function, alpha can
     % be determined by a simple formula.
    num = dot(df, s);
    den = (s.' * H * s);
    alpha = - num/den;
    %%%%%%%%%%%%%%%%%%%%
 elseif functionID == 2     
     desc_dir = 1;      % 1 if the given direcion is indeed a descent direction
     alpha_min = 0;
     h = 10;
     epsilon = 1e-8;
     rho = 0.5;
     if(dot(df,s)>0)
%          error('The given direction is not a descent direction');
         desc_dir = 0;
         s = - s;
     end
     GradF = getSens(x_init + (h)*s, functionID);
     PhiPrime = dot(GradF, s);  % directional derivative
     while(PhiPrime <= 0)       % Find initial boundary
         alpha_min = h;
         h = 2*h;  
         GradF = getSens(x_init + (h)*s, functionID);     
         PhiPrime = dot(GradF, s);
     end
     alpha_max = h;
     AlphaC = rho*alpha_min +  (1-rho)*alpha_max;
     GradF = getSens(x_init + AlphaC*s, functionID);
     PhiPrime = dot(GradF, s);
     while(abs(PhiPrime) > epsilon) % Bissection algorithm
         if(dot(getSens(x_init + alpha_max*s, functionID),s)*dot(getSens(x_init + alpha_min*s, functionID),s)>0)
             txt = sprintf('\n\n %s \n\n','Error in the bissection algorithm');
             error(txt);
         end
          if(n_step>100)
            warning('Precision of bissection was not reached for 1 iter');
            break;
         end
         n_step = n_step+1;
         AlphaC = rho*alpha_min +  (1-rho)*alpha_max;
         GradF = getSens(x_init + AlphaC*s, functionID);
         PhiPrime = dot(GradF, s);
         abs(PhiPrime);              
         if(PhiPrime < 0)
              alpha_min = AlphaC;
         elseif(PhiPrime > 0)
              alpha_max = AlphaC;
         elseif(PhiPrime==0)
             break;
         end
     end
     if(desc_dir==1)
        alpha = AlphaC;
     else
        alpha = -AlphaC;   
     end


     %%%%%%%%%%%%%%%%%%%%
     elseif functionID == 12
     functionID = 1
     desc_dir = 1;
     alpha_min = 0;
     h = 10;
     epsilon = 1e-8;
     rho = 0.5;
     if(dot(df,s)>0)
         desc_dir = 0;
         s = - s;
     end
     GradF = getSens(x_init + (h)*s, functionID);
     PhiPrime = dot(GradF, s);
     while(PhiPrime <= 0)
         alpha_min = h;
         h = 2*h;  
         GradF = getSens(x_init + (h)*s, functionID);     
         PhiPrime = dot(GradF, s);
     end
     alpha_max = h;
     AlphaC = rho*alpha_min +  (1-rho)*alpha_max;
     GradF = getSens(x_init + AlphaC*s, functionID);
     PhiPrime = dot(GradF, s);
     while(abs(PhiPrime) > epsilon)
         if(dot(getSens(x_init + alpha_max*s, functionID),s)*dot(getSens(x_init + alpha_min*s, functionID),s)>0)
             txt = sprintf('\n\n %s \n\n','Error in the bissection algorithm');
             error([txt]);
         end
         if(n_step>100 && abs(PhiPrime)< epsilon*1e3)
            warning('Precision of bissection was not reached for 1 iter');
            break;
         end
         n_step = n_step+1;
         AlphaC = rho*alpha_min +  (1-rho)*alpha_max;
         GradF = getSens(x_init + AlphaC*s, functionID);
         PhiPrime = dot(GradF, s);
         abs(PhiPrime);              
         if(PhiPrime < 0)
              alpha_min = AlphaC;
         elseif(PhiPrime > 0)
              alpha_max = AlphaC;
         end
     end
     if(desc_dir==1)
        alpha = AlphaC;
     else
        alpha = -AlphaC;   
     end
 end
 %%% you can use reshape to only consider vector columns
