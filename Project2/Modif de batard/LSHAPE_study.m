function LSHAPE_study(study) % 1 without projection mesh, 2 with projection %without projection vfrac
%% Base parameter
nx0 = 30;
ny0 = nx0;
volfrac = 0.5;
penal = 3.0;
rmin = 1.1;
ft = 1;
penals = [1 3 5 10];
rmin_study = [0.8 1.5 2 2.4];

%% Solving the problem
for i =1:4
    if(study==1)
    	[xPhys, Mnd, loop, Compliance] = top88(2^(i-1)*nx0,2^(i-1)*ny0,volfrac,penal,2^(i-1)*nx0*0.04,ft);
    elseif(study==2)
        [xPhys, Mnd, loop, Compliance] =  top110(2^(i-1)*nx0,2^(i-1)*ny0,volfrac,penal,2^(i-1)*nx0*0.05,3);
    elseif(study==3 )
        [xPhys, Mnd, loop, Compliance] = top88(4*nx0,4*ny0,volfrac-0.3+(i*0.1),penal,4*nx0*0.05,ft);  
    end
    xPhys_plot{i} = xPhys;
    Mnd_plot{i} = Mnd;
    loop_plot{i} = loop;
    Compliance_plot{i} = Compliance;
end

%% Plotting the results
Figure1=figure(1);clf;set(Figure1,'defaulttextinterpreter','latex');
hold on;
for i =1:4
subplot(2,2,i)
colormap(gray); imagesc(1-xPhys_plot{i}); caxis([0 1]); axis equal; axis off; drawnow;
end
end