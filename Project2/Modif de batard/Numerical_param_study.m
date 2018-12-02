function Numerical_param_study(study)  % study = 1 for mesh refinement; 2 for radius size; 3 for penalization factor; 4 for rmin+refinement
%% Base parameter
nx0 = 30;
ny0 = nx0/3;
volfrac = 0.5;
penal = 3.0;
rmin = 1.1;
ft = 1;
penals = [1 3 5 10];
rmin_study = [0.8 1.5 2 2.4];
%% Impact of the mesh refinement
for i =1:4
    if(study==1)
        [xPhys, Mnd, loop, Compliance] = top88(2^(i-1)*nx0,2^(i-1)*ny0,volfrac,penal,rmin,ft);
    elseif(study == 2)
        [xPhys, Mnd, loop, Compliance] = top88(2*nx0,2*ny0,volfrac,penal,rmin_study(i),ft);
    elseif(study == 3)
        ft = 2;
        [xPhys, Mnd, loop, Compliance] = top88(2*nx0,2*ny0,volfrac,penals(i),0.04*2*nx0,ft);  
    elseif(study==4)
        [xPhys, Mnd, loop, Compliance] = top88(2^(i-1)*nx0,2^(i-1)*ny0,volfrac,penal,0.04*2^(i-1)*nx0,ft);
    end
    xPhys_plot{i} = xPhys;
    Mnd_plot{i} = Mnd;
    loop_plot{i} = loop;
    Compliance_plot{i} = Compliance;
end

% Plot structure
Figure1=figure(1);clf;set(Figure1,'defaulttextinterpreter','latex');
hold on;
% set(gca,'fontsize',16,'fontname','Times','LineWidth',0.5);
% colormap(gray); imagesc(1-xPhys_plot{1}); caxis([0 1]); axis equal; axis off; 
% 
% Figure2=figure(2);clf;set(Figure2,'defaulttextinterpreter','latex');
% hold on;
% set(gca,'fontsize',16,'fontname','Times','LineWidth',0.5);
% colormap(gray); imagesc(1-xPhys_plot{2}); caxis([0 1]); axis equal; axis off; 
% 
% Figure3=figure(3);clf;set(Figure3,'defaulttextinterpreter','latex');
% hold on;
% set(gca,'fontsize',16,'fontname','Times','LineWidth',0.5);
% colormap(gray); imagesc(1-xPhys_plot{3}); caxis([0 1]); axis equal; axis off; 
% Figure4=figure(4);clf;set(Figure4,'defaulttextinterpreter','latex');
% hold on;
% set(gca,'fontsize',16,'fontname','Times','LineWidth',0.5);
% colormap(gray); imagesc(1-xPhys_plot{4}); caxis([0 1]); axis equal; axis off; 


for i =1:4
subplot(2,2,i)
colormap(gray); imagesc(1-xPhys_plot{i}); caxis([0 1]); axis equal; axis off; drawnow;
end

% Plot compliance
Figure2=figure(2);clf;set(Figure2,'defaulttextinterpreter','latex');
hold on;
set(gca,'fontsize',16,'fontname','Times','LineWidth',0.5);
for i =1:4
subplot(2,2,i)
plot(Compliance_plot{i},'r','linewidth',2);
axis([0 loop_plot{4} 0 max([max(Compliance_plot{1}) max(Compliance_plot{2}) max(Compliance_plot{3}) max(Compliance_plot{4})])]);
grid;
end

%% Mnd as a function of the refinement
if(study==1)
    i = 1;
    n = 12:12:120;
    nx0 = 12;
   while(nx0<=120)
       [xPhys, Mnd, loop, Compliance] = top88(nx0,nx0/3,volfrac,penal,rmin,ft);
       Mnd_plot_de_ouf(i) = Mnd;
       nx0 = nx0 +12;
       i = i+1;
   end
   Figure3=figure(3);clf;set(Figure3,'defaulttextinterpreter','latex');
   hold on;
   set(gca,'fontsize',25,'fontname','Times','LineWidth',0.5);
   plot(n,Mnd_plot_de_ouf,'r','linewidth',3);
   ylabel('Mnd');
   xlabel('Number of element in the x direction')
   grid;
end

%% Mnd as a function of the filter radius
if(study==2)
    rmin = 0.7
    i = 1;
    r = 0.7:0.05:2.5
   while(rmin<=2.5)
       [xPhys, Mnd, loop, Compliance] = top88(4*nx0,4*ny0,volfrac,penal,rmin,ft);
       Mnd_plot_de_ouf(i) = Mnd;
       rmin = rmin +0.05;
       i = i+1;
   end
   Figure3=figure(3);clf;set(Figure3,'defaulttextinterpreter','latex');
   hold on;
   set(gca,'fontsize',25,'fontname','Times','LineWidth',0.5);
   plot(r,Mnd_plot_de_ouf,'r','linewidth',3);
   ylabel('Mnd');
   xlabel('$r_{min}$')
   grid;
end

%% Mnd as a function of the penal value
if(study==3)
    penals = [1 2 3 4 5 8 9 10];
   for i = 1:5
       [xPhys, Mnd, loop, Compliance] = top88(2*nx0,2*ny0,volfrac,penals(i),0.04*2*nx0,ft);
       Mnd_plot_de_ouf(i) = Mnd;
   end
   for i = 1:8
       [xPhys, Mnd, loop, Compliance] = top88(2*nx0,2*ny0,volfrac,penals(i),0.04*2*nx0,ft);
       Mnd_plot_de_ouf(i) = Mnd;
   end
   Figure3=figure(3);clf;set(Figure3,'defaulttextinterpreter','latex');
   hold on;
   set(gca,'fontsize',16,'fontname','Times','LineWidth',0.5);
   plot(penals,Mnd_plot_de_ouf,'r','linewidth',2);
   ylabel('Mnd');
   xlabel('$Penality factor$')
   grid;
   
   [xPhys, Mnd, loop, Compliance] = top88(4*nx0,4*ny0,volfrac,1,0.04*4*nx0,ft);
   Compliance_plot_de_ouf{1} = Compliance;
   [xPhys, Mnd, loop, Compliance] = top88(4*nx0,4*ny0,volfrac,3,0.04*4*nx0,ft);
   Compliance_plot_de_ouf{2} = Compliance;
   
   
    Figure4=figure(4);clf;set(Figure4,'defaulttextinterpreter','latex');
    hold on;
    set(gca,'fontsize',60,'fontname','Times','LineWidth',0.5);
    subplot(2,1,1)
    hold on;    
    plot(Compliance_plot_de_ouf{1},'r','linewidth',2);  
    plot([30 30], [0 max([length(Compliance_plot_de_ouf{1}) length(Compliance_plot_de_ouf{2})])], '-.k','linewidth',2);
    plot([60 60], [0 max([length(Compliance_plot_de_ouf{1}) length(Compliance_plot_de_ouf{2})])], '-.k','linewidth',2)
    plot([90 90], [0 max([length(Compliance_plot_de_ouf{1}) length(Compliance_plot_de_ouf{2})])], '-.k','linewidth',2)
    plot([120 120], [0 max([length(Compliance_plot_de_ouf{1}) length(Compliance_plot_de_ouf{2})])], '-.k','linewidth',2)
    grid;
    ty = ylabel('Compliance');
    ty.FontSize = 20;
    tx  = xlabel('Iteration');
    tx.FontSize = 20;
    axis([0 max([length(Compliance_plot_de_ouf{1}) length(Compliance_plot_de_ouf{2})]) 0 max(Compliance_plot_de_ouf{1})]);
    set(gca,'FontSize',30)
    subplot(2,1,2)
    hold on;
    
    plot(Compliance_plot_de_ouf{2},'r','linewidth',2);
    grid;
    ty = ylabel('Compliance');
    ty.FontSize = 20;
    tx  = xlabel('Iteration');
    tx.FontSize = 20;
    set(gca,'FontSize',30)
    axis([0 max([length(Compliance_plot_de_ouf{1}) length(Compliance_plot_de_ouf{2})]) 0 max(Compliance_plot_de_ouf{2})]);
end

end