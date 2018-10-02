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


function [OptSol,FVal] = TD1_2018(xinit)

disp('Which optimization method do you want to use? Press:')
disp('     1 for Steepest descent method')
disp('     21 for Conjugate gradient method')
disp('     22 for Conjugate direction method using Fletcher-Reeves approximation')
disp('     3 for Newton method')
disp('     4 for Quasi-Newton method - Boyden-Fletcher-Goldfard-Shanno (BFGS)')
prompt = '';
ind = input(prompt);

% Parameters --------------------------------------------------------------------------
functionID = 1;             % =1 if min f(x,y)= 2x^2 + 2y^2 - 3xy + 2y^2 - 2x + 10y - 1
                            % =2 if min f(x,y)= 2x^4 + 2y^2 - 3xy + 2y^2 - 2x + 10y - 1
MaxIter=25;                 %Maximum number of iterations
n=2;                        %Dimension of the problem
xinit=reshape(xinit,2,1);   %Be sure it is a column vector
Epsilon=1e-4;               %Tolerance
x=zeros(n,MaxIter);         %Initialization of x
x(:,1)=xinit;               %Put xinit in vector x.

% Methods ------------------------------------------------
if ind==1
    disp(' You chose the steepest descent method.');
    
    for i=1:MaxIter
           H=getHess(x(:,i), functionID);
        df = getSens(x(:,i),functionID);
        s = -df;
        Norm_df = norm(df);
        alpha = getalpha(x(:,i),s,df,H,functionID);
        x(:,i+1) = x(:,i) +(alpha*s)/Norm_df
    end
    x=x(:,1:i); %Remove the zero elements due to the initialization step
    
elseif ind==21
    disp(' You chose the conjugate gradient method.')
    H=getHess();
    
    for i=1:MaxIter
       %%%% ADD YOUR CODE
    end
    x=x(:,1:i); %Remove the zero elements due to the initialization step

elseif ind==22
    disp(' You chose the conjugate direction method using Fletcher-Reeves approximation.')
    H=getHess();
    
    for i=1:MaxIter
       %%%% ADD YOUR CODE
    end
    x=x(:,1:i); %Remove the zero elements due to the initialization step

elseif ind==3
    disp(' You chose the Newton method.')
    H=getHess();
    
    for i=1:MaxIter
        %%%% ADD YOUR CODE
    end
    x=x(:,1:i); %Remove the zero elements due to the initialization step
    	
elseif ind==4
    disp(' You chose the quasi-Newton method (Boyden-Fletcher-Goldfard-Shanno (BFGS)).')
    Sinit=eye(n);
    H=getHess(); %For the exact line search
    
    for i=1:MaxIter
       %%%% ADD YOUR CODE
    end
    x=x(:,1:i); %Remove the zero elements due to the initialization step
     
else
    disp('Error: Your choice does not match the proposed methods.')
    return
end

OptSol=x(:,end);
disp(['The optimal point is: x = ' num2str(OptSol(1)) ', y = ' num2str(OptSol(2)) '.'])
FVal=getObjFVal(OptSol,functionID);
disp(['The objective function value is: ' num2str(FVal) '.'])

%Plot the function with the optimization path.
lb=-3;
up=3;
xi=lb:0.1:up;
f=zeros(length(xi),length(xi));
for i=1:length(xi)
    for j=1:length(xi)
        f(j,i)=getObjFVal([xi(i);xi(j)],functionID);
    end
end
figure('name','Optimization path')
hold on
axis equal
xlabel('x_1')
ylabel('x_2')
axis([lb up lb up])
title('Optimization path')
[C,h]=contour(xi,xi,f,[-4:2:8 10:10:50 75:25:200]);
clabel(C,h);
hold on
ind=1;
for i=1:size(x,2)-1
    plot(x(1,i),x(2,i),'.c','markersize',30)
    plot([x(1,i) x(1,i+1)],[x(2,i) x(2,i+1)],'c','linewidth',2)
    text(x(1,i),x(2,i),num2str(ind-1),'horizontalalignment','center','verticalalignment','middle')
    ind=ind+1;
end
plot(x(1,end),x(2,end),'.c','markersize',30)
text(x(1,end),x(2,end),num2str(ind-1),'horizontalalignment','center','verticalalignment','middle')


function fval = getObjFVal(x,functionID)
 if functionID == 1
     fval = 2*(x(1))^2-3*(x(1)*x(2))+2*(x(2))^2 - 2*x(1)+10*x(2)-1;
 elseif functionID == 2
     fval = 2*x(1)^4 - 3*x(1)*x(2) + 2*x(2)^2 - 2*x(1) + 10*x(2) - 1;
 end

 function df = getSens(x,functionID)
 if functionID == 1
     df(1) = 4*x(1) - 3*x(2) - 2;
     df(2) = -3*x(1) + 4*x(2) + 10
     df = df.';
 elseif functionID == 2
     %%%% ADD YOUR CODE
 end

 function H = getHess(x,functionID)
 if functionID == 1
     H(1,1) = 4;
     H(1,2) = -3;
     H(2,1) = -3;
     H(2,2) = 4;
 elseif functionID == 2
     %%%% ADD YOUR CODE
 end

 function alpha = getalpha(x_init,s,df,H,functionID)
 if functionID == 1
     %df = getSens(x_init);
     %H = getHess(x_init);
    num = dot(df, s);
    den = (s.' * H * s);
    alpha = - num/den;
 elseif functionID == 2
     %%%% ADD YOUR CODE
 end
 %%% you can use reshape to only consider vector columns
