clear all
close all
nelx = 240;
nely = round(nelx/2);
volfrac = 0.5;
penal = 3;
ft = 1;
rmin = 3.2;
mkdir('TwoLoad4');
disp('Running the two load case study');
%% CASE 1:
WhichLoading = 6;
WhichBoundary = 3;
if (mod(nelx, 3) ~= 0)
    disp('hello')
end
tic
[xPhysAlpha1, MndAlpha1, loopAlpha1, ComplianceAlpha1, SvmAlpha1] = top88BoundaryLoading(nelx,nely,volfrac,penal,rmin,ft, WhichLoading, WhichBoundary);
disp('done case 1, plotting...');
toc
IterAlpha1 = 1:loopAlpha1;
cd('TwoLoad4');
figure(2);
CompPlot = myPlot(IterAlpha1, ComplianceAlpha1, 'Iteration', 'Compliance');
matlab2tikz('ComplianceCase1.tex','width', '0.8\textwidth', 'height', '0.4\textwidth');
figure(3);
colormap(gray); imagesc(1-xPhysAlpha1); caxis([0 1]); axis equal; axis off;
%hold on; myPlot(0,0);
print('OptimisedGeomCase1','-depsc');
cd('..');

%% CASE 2:
WhichLoading = 7;
WhichBoundary = 3;
if (mod(nelx, 3) ~= 0)
    disp('hello')
end
tic
[xPhysAlpha1, MndAlpha1, loopAlpha1, ComplianceAlpha1, SvmAlpha1] = top88BoundaryLoading(nelx,nely,volfrac,penal,rmin,ft, WhichLoading, WhichBoundary);
disp('done case 1, plotting...');
toc
IterAlpha1 = 1:loopAlpha1;
cd('TwoLoad4');
figure(2);
CompPlot = myPlot(IterAlpha1, ComplianceAlpha1, 'Iteration', 'Compliance');
matlab2tikz('ComplianceCase2.tex','width', '0.8\textwidth', 'height', '0.4\textwidth');
figure(3);
colormap(gray); imagesc(1-xPhysAlpha1); caxis([0 1]); axis equal; axis off;
%hold on; myPlot(0,0);
print('OptimisedGeomCase2','-depsc');
cd('..');

% %% CASE 3:
% WhichLoading = 8;
% WhichBoundary = 3;
% if (mod(nelx, 3) ~= 0)
%     disp('hello')
% end
% tic
% [xPhysAlpha3, MndAlpha3, loopAlpha3, ComplianceAlpha3, SvmAlpha3] = top88BoundaryLoading(nelx,nely,volfrac,penal,rmin,ft, WhichLoading, WhichBoundary);
% disp('done case 3, plotting...');
% toc
% IterAlpha3 = 1:loopAlpha3;
% cd('TwoLoad3');
% figure(2);
% CompPlot = myPlot(IterAlpha3, ComplianceAlpha3, 'Iteration', 'Compliance');
% matlab2tikz('ComplianceCase3.tex','width', '0.8\textwidth', 'height', '0.4\textwidth');
% figure(3);
% colormap(gray); imagesc(1-xPhysAlpha3); caxis([0 1]); axis equal; axis off;
% %hold on; myPlot(0,0);
% print('OptimisedGeomCase3','-depsc');
% cd('..');
% %% CASE 3:
% WhichLoading = 9;
% WhichBoundary = 3;
% if (mod(nelx, 3) ~= 0)
%     disp('hello')
% end
% tic
% [xPhysAlpha33, MndAlpha33, loopAlpha33, ComplianceAlpha33, SvmAlpha33] = top88BoundaryLoading(nelx,nely,volfrac,penal,rmin,ft, WhichLoading, WhichBoundary);
% disp('done case 3, plotting...');
% toc
% IterAlpha33 = 1:loopAlpha33;
% cd('TwoLoadbis');
% figure(2);
% CompPlot = myPlot(IterAlpha33, ComplianceAlpha33, 'Iteration', 'Compliance');
% matlab2tikz('ComplianceCase33.tex','width', '0.8\textwidth', 'height', '0.4\textwidth');
% figure(3);
% colormap(gray); imagesc(1-xPhysAlpha33); caxis([0 1]); axis equal; axis off;
% %hold on; myPlot(0,0);
% print('OptimisedGeomCase33','-depsc');
% cd('..');