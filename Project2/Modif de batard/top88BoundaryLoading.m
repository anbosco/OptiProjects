%%%% AN 88 LINE TOPOLOGY OPTIMIZATION CODE Nov, 2010 %%%%
function [xPhys, Mnd, loop, Compliance, Svm] = top88BoundaryLoading(nelx,nely,volfrac,penal,rmin,ft, WhichLoading, WhichBoundary)
close all;
Cont =1;
dp = 0.5;
nelx
if (nelx/nely ~= 2 || mod(nelx, 3) ~= 0 )
    error('This case is defined for nelx/nely = 2 and mod(nelx, 3) == 0');
end
%% MATERIAL PROPERTIES
E0 = 1;
Emin = 1e-9;
nu = 0.3;
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
%% L-SHAPE
% METTRE LSHAPE = 1 pour activer le cas avec LSHAPE
% J'ai fait �a parce que sinon fallait commenter et d�comenter
% � plusieurs endroit dans le code (cfr. vers ligne 138)
LSHAPE = 0;
if LSHAPE ==0
    % DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
    %Load
    if WhichLoading == 1
        F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1); % CAS DE BASE
        U = zeros(2*(nely+1)*(nelx+1),1);
    elseif WhichLoading == 2
        F = sparse(1+(nely + 1)*(nelx + 1) - nely ,1,-2, 2*(nely+1)*(nelx+1),1);
        U = zeros(2*(nely+1)*(nelx+1),1);         % CAS FORCE EQUIVALENT
    elseif WhichLoading == 3
        F = sparse(round((2/3)*(nely+1)*(nelx) + 2) ,1,1, 2*(nely+1)*(nelx+1),1);
        F(round((4/3)*(nely+1)*(nelx) + 2) ,1) = -1;
        U = zeros(2*(nely+1)*(nelx+1),1);           % CAS 2 FORCES EN MEME TEMPS
     elseif WhichLoading == 33
        F = sparse(round((2/3)*(nely+1)*(nelx) + 2) ,1,1, 2*(nely+1)*(nelx+1),1);
        F(round((4/3)*(nely+1)*(nelx) + 2) ,1) = 1;
        U = zeros(2*(nely+1)*(nelx+1),1);    
    elseif WhichLoading == 4 % Increasing magnitude
        F = sparse(round((2/3)*(nely+1)*(nelx) + 2) ,1,10, 2*(nely+1)*(nelx+1),1);
        F(round((4/3)*(nely+1)*(nelx) + 2) ,1) = 10;
        U = zeros(2*(nely+1)*(nelx+1),1);
    elseif WhichLoading == 5 % Direction;
        F = sparse(round((2/3)*(nely+1)*(nelx) + 2) ,1,-1, 2*(nely+1)*(nelx+1),1);
        F(round((4/3)*(nely+1)*(nelx) + 2) ,1) = -1;
        U = zeros(2*(nely+1)*(nelx+1),1);
    elseif WhichLoading == 6
        F = sparse(2*(nely+1)*(nelx+1),2);
        F(round((2/3)*(nely+1)*(nelx) + 2),1) = 1;
        F(round((4/3)*(nely+1)*(nelx) + 2),2) = 1;
        U = zeros(2*(nely+1)*(nelx+1),2);         % CAS 2 FORCES PAS EN MEME TEMPS
     elseif WhichLoading == 7
        F = sparse(2*(nely+1)*(nelx+1),2);
        F(round((2/3)*(nely+1)*(nelx) + 2),1) = 1;
        F(round((4/3)*(nely+1)*(nelx) + 2),2) = -1;
        U = zeros(2*(nely+1)*(nelx+1),2);
      elseif WhichLoading == 8
        F = sparse(2*(nely+1)*(nelx+1),2);
        F(round((2/3)*(nely+1)*(nelx) + 2),1) = 1;
        F(round((4/3)*(nely+1)*(nelx) + 2),2) = 0;
        U = zeros(2*(nely+1)*(nelx+1),2);
              elseif WhichLoading == 9
        F = sparse(2*(nely+1)*(nelx+1),2);
        F(round((2/3)*(nely+1)*(nelx) + 2),1) = 10;
        F(round((4/3)*(nely+1)*(nelx) + 2),2) = 10;
        U = zeros(2*(nely+1)*(nelx+1),2);
    else
        error('Use a valid loading ID');
    end
    
    % Craquage sur les loading because i can :D 
    
    %BC
    if WhichBoundary == 1
        fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
        alldofs = [1:2*(nely+1)*(nelx+1)];
        freedofs = setdiff(alldofs,fixeddofs);      % CAS DE BASE
    elseif WhichBoundary == 2
        fixeddofs = union([2*(nely+1)-1:2*(nely+1)],[2*(nelx+1)*(nely+1)-1:2*(nelx+1)*(nely+1)]);
        alldofs = [1:2*(nely+1)*(nelx+1)];
        freedofs = setdiff(alldofs,fixeddofs);      % CAS 2
    elseif WhichBoundary == 3
        fixeddofs = union([2*(nely+1)-1:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
        alldofs = [1:2*(nely+1)*(nelx+1)];
        freedofs = setdiff(alldofs,fixeddofs);      % CAS 3
    else
        error('Use a valid boundary ID');
    end
end

%% L-SHAPE
% METTRE LSHAPE = 1 pour activer le cas avec LSHAPE
% J'ai fait �a parce que sinon fallait commenter et d�commenter
% � plusieurs endroit dans le code (cfr. vers ligne 138)
if LSHAPE == 1
    passive = zeros(nely,nelx);
    for i = 1:nelx
        for j = 1:nely
            if (i > round(0.4*nelx) && j < (0.6*nely))
                passive(j,i) = 1;
            end
        end
    end
    
    % BC L-SHAPE
    % A MODIFIER, ELLE EST FAUSSE :(
    % fixeddofs = union([[0:2:(0.4*nelx+2)].*(nelx+1) + 2 + 2*(nely+1)],[2]);
    fixeddofs = union([2:2*nely+2:0.4*nelx*2*nely+1+2*nely+2],[1:2*nely+2:0.4*nelx*2*nely+1+2*nely+1]);
    
    alldofs = [1:2*(nely+1)*(nelx+1)];
    freedofs = setdiff(alldofs,fixeddofs);
    
    F = sparse(2*(nelx+1)*(nely+1)-0.4*nely*2 ,1,1/(0.1*nely), 2*(nely+1)*(nelx+1),1);
    for i = 2*(nelx+1)*(nely+1)-0.4*nely*2+2:2:2*(nelx+1)*(nely+1)-0.3*nely*2
        F(i,1) = 1/(0.1*nely);
    end
    U = zeros(2*(nely+1)*(nelx+1),1);           % CAS 2 FORCES EN MEME TEMPS
end

% Plot BC
%gfix(nelx,nely,fixeddofs,F,[])
%figure;
%  error('On fait les BC putain')
%error('On fait les BC putain')

%% PREPARE FILTER
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
    for j1 = 1:nely
        e1 = (i1-1)*nely+j1;
        for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
            for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                e2 = (i2-1)*nely+j2;
                k = k+1;
                iH(k) = e1;
                jH(k) = e2;
                sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
            end
        end
    end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
%% INITIALIZE ITERATION
x = repmat(volfrac,nely,nelx);
xPhys = x;
loop = 0;
change = 1;
%% START ITERATION
while change > 0.01
    loop = loop + 1;
    %% FE-ANALYSIS
    sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
    K = sparse(iK,jK,sK); K = (K+K')/2;
    U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    c = 0;
    dc = 0;
    for(i=1:size(F,2))
        Ui = U(:,i);
        ce = reshape(sum((Ui(edofMat)*KE).*Ui(edofMat),2),nely,nelx);
        c = c + sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));
        dc = dc - penal*(E0-Emin)*xPhys.^(penal-1).*ce;
    end
    dv = ones(nely,nelx);
    %% FILTERING/MODIFICATION OF SENSITIVITIES
    if ft == 1
        dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
    elseif ft == 2
        dc(:) = H*(dc(:)./Hs);
        dv(:) = H*(dv(:)./Hs);
    end
    %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
    l1 = 0; l2 = 1e9; move = 0.2;
    %%
    while (l2-l1)/(l1+l2) > 1e-3
        lmid = 0.5*(l2+l1);
        xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));
        
        if ft == 1
            xPhys = xnew;
        elseif ft == 2
            xPhys(:) = (H*xnew(:))./Hs;
        end
        if LSHAPE == 1
            xPhys(passive==1) = 0;
            xPhys(passive==2) = 1;
        end
        if sum(xPhys(:)) > volfrac*nelx*nely, l1 = lmid; else l2 = lmid; end
    end
    change = max(abs(xnew(:)-x(:)));
    x = xnew;
    
    % Coucou je suis pas efficace et je nique ton code
    temp = 4*x.*(1-x);
    Mnd = sum(sum(temp))/(length(x(:,1))*length(x(1,:)));
    Mnd = Mnd*100;
    %% PRINT RESULTS
    %   fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c, ...
    %     mean(xPhys(:)),change);
    % PLOT DENSITIES
    %   colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
    %% Continuation
    if(Cont==1 && penal< 3 && mod(loop,30)==0)
        penal = penal+dp;
    end
    %% Saving compliance
    Compliance(loop)=c;
end
[Svm] = computestress(U,edofMat,E0,nu,penal,xPhys);
% gfix(nelx,nely,fixeddofs,F,[])
% figure;
% plot(Compliance);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by E. Andreassen, A. Clausen, M. Schevenels,%
% B. S. Lazarov and O. Sigmund,  Department of Solid  Mechanics,           %
%  Technical University of Denmark,                                        %
%  DK-2800 Lyngby, Denmark.                                                %
% Please sent your comments to: sigmund@fam.dtu.dk                         %
%                                                                          %
% The code is intended for educational purposes and theoretical details    %
% are discussed in the paper                                               %
% "Efficient topology optimization in MATLAB using 88 lines of code,       %
% E. Andreassen, A. Clausen, M. Schevenels,                                %
% B. S. Lazarov and O. Sigmund, Struct Multidisc Optim, 2010               %
% This version is based on earlier 99-line code                            %
% by Ole Sigmund (2001), Structural and Multidisciplinary Optimization,    %
% Vol 21, pp. 120--127.                                                    %
%                                                                          %
% The code as well as a postscript version of the paper can be             %
% downloaded from the web-site: http://www.topopt.dtu.dk                   %
%                                                                          %
% Disclaimer:                                                              %
% The authors reserves all rights but do not guaranty that the code is     %
% free from errors. Furthermore, we shall not be liable in any event       %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

