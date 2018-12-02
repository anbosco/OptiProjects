function LSHAPE_study
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
    [xPhys, Mnd, loop, Compliance] = top88(2^(i-1)*nx0,2^(i-1)*ny0,volfrac,2^(i-1)*nx0*0.04,rmin,ft);
    xPhys_plot{i} = xPhys;
    Mnd_plot{i} = Mnd;
    loop_plot{i} = loop;
    Compliance_plot{i} = Compliance;
end
Figure1=figure(1);clf;set(Figure1,'defaulttextinterpreter','latex');
hold on;
for i =1:4
subplot(2,2,i)
colormap(gray); imagesc(1-xPhys_plot{i}); caxis([0 1]); axis equal; axis off; drawnow;
end
end