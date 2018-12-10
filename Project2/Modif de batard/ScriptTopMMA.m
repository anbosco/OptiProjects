clear all;
close all;
nelx = 120;
nely = round(nelx/2);
volfrac = 0.5;
penal = 3;
ft = 1;
rmin = 1.5;
WhichBoundary = 3;
WhichLoading = 7;
%% Case MMA
tic
[xPhys, Mnd, loop, Compliance] = topMMA(nelx,nely,volfrac,penal,rmin,ft);
time{1} = toc
    xPhys_plot{1} = xPhys;
    Mnd_plot{1} = Mnd;
loop_plot{1} = loop;
    Compliance_plot{1} = Compliance;
%% Case top88
tic
    [xPhys, Mnd, loop, Compliance, Svm] = ...
        top88BoundaryLoading(nelx,nely,volfrac,penal,rmin,ft, WhichLoading,...
        WhichBoundary);
time{2} = toc    
        xPhys_plot{2} = xPhys;
    Mnd_plot{2} = Mnd;
    loop_plot{2} = loop;
    Compliance_plot{2} = Compliance;
    %% Stuff
    mkdir('TopMMAbis');
    cd('TopMMAbis');
%% Case 1:
Iter1 = 1:loop_plot{1};
figure(1);
CompPlot = myPlot(Iter1, Compliance_plot{1}, 'Iteration', 'Compliance');
matlab2tikz('ComplianceCase1.tex','width', '0.8\textwidth', 'height', '0.4\textwidth');
figure(11);
colormap(gray); imagesc(1-xPhys_plot{1}); caxis([0 1]); axis equal;axis off;
%hold on; myPlot(0,0);
print('OptimisedGeomCase1','-depsc');

%% Case 2:
Iter2 = 1:loop_plot{2};
figure(2);
CompPlot = myPlot(Iter2, Compliance_plot{2}, 'Iteration', 'Compliance');
matlab2tikz('ComplianceCase2.tex','width', '0.75\textwidth', 'height', '0.35\textwidth');
figure(22);
colormap(gray); imagesc(1-xPhys_plot{2}); caxis([0 1]);axis equal ;axis off;
%hold on; myPlot(0,0);
print('OptimisedGeomCase2','-depsc');

%%
figure;
CompPlot = myPlot(Iter1, Compliance_plot{1}, 'Iteration', 'Compliance');
hold on; CompPlot = myPlot(Iter2, Compliance_plot{2}, 'Iteration', 'Compliance');
legend('MMA', '$\texttt{top88}$', 'location', 'northeast');
matlab2tikz('ComplianceCase12.tex','width', '0.75\textwidth', 'height', '0.35\textwidth');
