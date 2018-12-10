%%%% AN 88 LINE TOPOLOGY OPTIMIZATION CODE Nov, 2010 %%%%
function [xPhys, Mnd, loop, Compliance] = topMMA(nelx,nely,volfrac,penal,rmin,ft)
close all;
Cont =1;
dp = 0.5;
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
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
%Load
% F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1); % CAS DE BASE
% U = zeros(2*(nely+1)*(nelx+1),1);

% F = sparse(1+(nely + 1)*(nelx + 1) - nely ,1,-2, 2*(nely+1)*(nelx+1),1);
% U = zeros(2*(nely+1)*(nelx+1),1);         % CAS FORCE EQUIVALENT

% F = sparse(round((2/3)*(nely+1)*(nelx) + 2) ,1,1, 2*(nely+1)*(nelx+1),1);
% F(round((4/3)*(nely+1)*(nelx) + 2) ,1) = 1;
% U = zeros(2*(nely+1)*(nelx+1),1);           % CAS 2 FORCES EN MEME TEMPS

F = sparse(2*(nely+1)*(nelx+1),2);
F(round((2/3)*(nely+1)*(nelx) + 2),1) = 1;
F(round((4/3)*(nely+1)*(nelx) + 2),2) = -1;
U = zeros(2*(nely+1)*(nelx+1),2);         % CAS 2 FORCES PAS EN MEME TEMPS

% BC
% fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
% alldofs = [1:2*(nely+1)*(nelx+1)];
% freedofs = setdiff(alldofs,fixeddofs);      % CAS DE BASE

% fixeddofs = union([2*(nely+1)-1:2*(nely+1)],[2*(nelx+1)*(nely+1)-1:2*(nelx+1)*(nely+1)]);
% alldofs = [1:2*(nely+1)*(nelx+1)];
% freedofs = setdiff(alldofs,fixeddofs);      % CAS 2

fixeddofs = union([2*(nely+1)-1:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);      % CAS 3

%% L-SHAPE  
% METTRE LSHAPE = 1 pour activer le cas avec LSHAPE 
% J'ai fait �a parce que sinon fallait commenter et d�comenter
% � plusieurs endroit dans le code (cfr. vers ligne 138)

% Plot BC
gfix(nelx,nely,fixeddofs,F,[])
figure;
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
V0 = sum(sum(x));
xPhys = x;
loop = 0;
change = 1;
%% Init MMA
  nVar = nelx*nely;
  xmin = zeros(nelx*nely, 1);
  xmax = ones(nelx*nely, 1);
  low = zeros(nelx*nely,1);
  upp = zeros(nelx*nely,1);
  xold1 = zeros(nelx*nely,1);
  xold2 = zeros(nelx*nely,1);
   df0dx2 = zeros(nelx*nely,1);
    dfdx2 = zeros(nelx*nely,1);
      a0mma = 1; amma = 0; cmma = 1000; dmma = 0;
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
%   while (l2-l1)/(l1+l2) > 1e-3
%     lmid = 0.5*(l2+l1);
%     xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));
%     iter = iter+1;    
%     if ft == 1
%       xPhys = xnew;
%     elseif ft == 2
%       xPhys(:) = (H*xnew(:))./Hs;
%     end
%     if LSHAPE == 1
%     xPhys(passive==1) = 0;
%     xPhys(passive==2) = 1;
%     end
%     if sum(xPhys(:)) > volfrac*nelx*nely, l1 = lmid; else l2 = lmid; end
%     break;
%   end
 
  %% MMA
 % fval = mean(xPhys(:)) - volfrac;
  fval = (sum(sum(xPhys)))/(nVar*volfrac) - 1;
  dfdx = (dv(:))/(nVar*volfrac);
  nConstr = 1;
  [xmma,~,~,~,~,~,~,~,~,low,upp] = ...
mmasub(nConstr,nVar,loop,x(:),xmin,xmax,xold1,xold2, ...
c,dc(:),df0dx2,fval,dfdx,dfdx2,low,upp,a0mma,amma,cmma,dmma);
    xold2 = xold1;
    xold1 = x(:);
    
    xnew = reshape(xmma, nely, nelx);
  
  change = max(abs(xnew(:)-x(:)));
  if ft == 1
      xPhys = xnew;
  elseif ft == 2
      xPhys(:) = (H*xnew(:))./Hs;
  end
  x = xnew;
%   % Coucou je suis pas efficace et je nique ton code
     temp = 4*x.*(1-x);
     Mnd = sum(sum(temp))/(length(x(:,1))*length(x(1,:)));
     Mnd = Mnd*100;
  %% PRINT RESULTS
 % fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c, ...
 %   mean(xPhys(:)),change);
  %% PLOT DENSITIES
 % colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
  %% Continuation
  if(Cont==1 && penal< 3 && mod(loop,30)==0)
     penal = penal+dp;      
  end
  %% Saving compliance
  Compliance(loop)=c;
end
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

