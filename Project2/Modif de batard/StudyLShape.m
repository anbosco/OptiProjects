clear all
close all
nx0 = 30;
ny0 = nx0;
volfrac = 0.5;
penal = 3.0;
rmin = 1.1;
ft = 1;
penals = [1 3 5 10];
rmin_study = [0.8 1.5 2 2.4];
study = 3;
%% Solving the problem
for i =1:4
    i
    if(study==1)
    	[xPhys, Mnd, loop, Compliance] = top88(2^(i-1)*nx0,2^(i-1)*ny0,volfrac,penal,2^(i-1)*nx0*0.04,ft);
    elseif(study==2)
        [xPhys, Mnd, loop, Compliance] =  top110(2^(i-1)*nx0,2^(i-1)*ny0,volfrac,penal,2^(i-1)*nx0*0.05,3);
    elseif(study==3 )
        [xPhys, Mnd, loop, Compliance, Svm] = top88(4*nx0,4*ny0,volfrac-0.3+(i*0.1),penal,4*nx0*0.05,ft);  
    end
    xPhys_plot{i} = xPhys;
    Mnd_plot{i} = Mnd;
    loop_plot{i} = loop;
    Compliance_plot{i} = Compliance;
    Svm_plot{i} = Svm;
end

%% Plotting the results
mkdir('LShape');
cd('LShape');
%% Case 1:
Iter1 = 1:loop_plot{1};
figure(1);
CompPlot = myPlot(Iter1, Compliance_plot{1}, 'Iteration', 'Compliance');
matlab2tikz('ComplianceCase1.tex','width', '0.8\textwidth', 'height', '0.4\textwidth');
figure(11);
colormap(gray); imagesc(1-xPhys_plot{1}); caxis([0 1]); axis square;axis off;
%hold on; myPlot(0,0);
print('OptimisedGeomCase1','-depsc');
figure(111);
colormap(jet); imagesc(Svm_plot{1}); caxis([0 1]);axis square; cb2=colorbar;cb2.TickLabelInterpreter = 'latex'; axis off;
%hold on; myPlot(0,0);
set(gca,'fontsize', 20);
print('SigmaVMCase1','-depsc');

%% Case 2:
Iter2 = 1:loop_plot{2};
figure(2);
CompPlot = myPlot(Iter2, Compliance_plot{2}, 'Iteration', 'Compliance');
matlab2tikz('ComplianceCase2.tex','width', '0.75\textwidth', 'height', '0.35\textwidth');
figure(22);
colormap(gray); imagesc(1-xPhys_plot{2}); caxis([0 1]);axis equal ;axis off;
%hold on; myPlot(0,0);
print('OptimisedGeomCase2','-depsc');
figure(222)
colormap(jet); imagesc(Svm_plot{2}); caxis([0 1]); axis square;cb2=colorbar;cb2.TickLabelInterpreter = 'latex'; axis off;
%hold on; myPlot(0,0);
set(gca,'fontsize', 20);
print('SigmaVMCase2','-depsc');

%% Case 3:
Iter3 = 1:loop_plot{3};
figure(3);
CompPlot = myPlot(Iter3, Compliance_plot{3}, 'Iteration', 'Compliance');
matlab2tikz('ComplianceCase3.tex','width', '0.75\textwidth', 'height', '0.35\textwidth');
figure(33);
colormap(gray); imagesc(1-xPhys_plot{3}); caxis([0 1]); axis square;axis off;
%hold on; myPlot(0,0);
print('OptimisedGeomCase3','-depsc');
figure(333)
colormap(jet); imagesc(Svm_plot{3}); caxis([0 1]); axis square;cb2=colorbar;cb2.TickLabelInterpreter = 'latex'; axis off;
%hold on; myPlot(0,0);
set(gca,'fontsize', 20);
print('SigmaVMCase3','-depsc');

%% Case 4:
Iter4 = 1:loop_plot{4};
figure(4);
CompPlot = myPlot(Iter4, Compliance_plot{4}, 'Iteration', 'Compliance');
matlab2tikz('ComplianceCase4.tex','width', '0.75\textwidth', 'height', '0.35\textwidth');
figure(44);
colormap(gray); imagesc(1-xPhys_plot{4}); caxis([0 1]); axis square;axis off;
%hold on; myPlot(0,0);
print('OptimisedGeomCase4','-depsc');
figure(444)
colormap(jet); imagesc(Svm_plot{4}); caxis([0 1]); axis square;cb2=colorbar;cb2.TickLabelInterpreter = 'latex'; axis off;
%hold on; myPlot(0,0);
set(gca,'fontsize', 20);
print('SigmaVMCase4','-depsc');
cd('..');