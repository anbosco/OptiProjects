function [] = gfix(nelx,nely,fixdof,FORCE,PE)
[nodo] = enum_nodos(nelx,nely);
c=0;
% ------------------ fixeddof ------------------------------------
for i=1:size(fixdof,2)
    c=c+1;
    [row,col]=find(nodo(:,(3:4))==fixdof(1,i));
    
    if col==1
        dire_fixed(c,:)=0.5*[1 0];
    else
        dire_fixed(c,:)=0.5*[0 1];
    end
    fixed_num(c,:)=nodo(row,(1:2));
end
%------------------------------------------------------------------

FORCE=full(FORCE);
ff=0.8;
for j=1:size(FORCE,2)
    c=0;
    for i=1:size(FORCE,1)
        if abs(FORCE(i,j))>0
            c=c+1;
            if mod(i,2)==1
                dire2{j,1}(c,:)=FORCE(i,j)*[1 0]*ff;
            else
                dire2{j,1}(c,:)=FORCE(i,j)*[0 1]*ff;
            end
            [row,col]=find(nodo(:,(3:4))==i);
            BC2{j,1}(c,:)=nodo(row,(1:2));
        end
    end
end


%---------------- plot creation ----------------------------------
hold off
figure(1)
plot(nodo(:,1),nodo(:,2),'k.','markers',3)
axis equal
axis tight
hold on
x1=nodo(1,1)                            ;   y1=nodo(1,2);
x2=nodo(nelx+1,1)                       ;   y2=nodo(nelx+1,2);
x3=nodo((nelx+1)*(nely+1),1)            ;   y3=nodo((nelx+1)*(nely+1),2);
x4=nodo((nelx+1)*(nely+1)-(nelx),1)     ;   y4=nodo((nelx+1)*(nely+1)-nelx,2);
plot([x1,x2,x3,x4,x1],[y1,y2,y3,y4,y1],'k--');
axis([-4 nelx+4 -4 nely+4]);
%---------------- plot fixed dof --------------------------------
if size(fixdof,2)>0
    
    for i=1:size(fixed_num,1)
        clear X1 X2 Y1 Y2
        X1=fixed_num(i,1);
        Y1=fixed_num(i,2);
        
        if X1<nelx
            X2=X1-dire_fixed(i,1);
        else
            X2=X1+dire_fixed(i,1);
        end
        
        if Y1<nely
            Y2=Y1-dire_fixed(i,2);
        else
            Y2=Y1-dire_fixed(i,2);
        end
        plot([X1 X2],[Y1 Y2],'b-','linewidth',2);
        plot(X2,Y2,'b.','markers',20);
    end
end
%---------------- plot loads --------------------------------
CO{1,1}='red'; CO{2,1}='green';  CO{3,1}='cyan'; CO{2,1}='magenta';
for j=1:size(FORCE,2)
    for i=1:size(BC2{j,1},1)
        lx=dire2{j,1}(i,1);
        ly=dire2{j,1}(i,2);
        if BC2{j,1}(i,2)==0
            x1=BC2{j,1}(i,1);
            y1=BC2{j,1}(i,2);
        else
            x1=BC2{j,1}(i,1)-lx;
            y1=BC2{j,1}(i,2)-ly;
        end

        p=quiver(x1,y1,lx,ly,0);
        set(p,'Color',strcat(CO{j,1}));
        set(p,'LineWidth',2.3);
        set(p,'MaxHeadSize',3);
    end
end
%---------------- plot passives --------------------------------
if sum(sum(PE))==0
else
    map = zeros(nely,nelx);
    map(PE)=1;
    for i = 1:nelx
        for j = 1:nely
            if map(j,i)~=0
                plot(i-0.5,j-0.5,'sr','MarkerSize',4);
            end
        end
    end
end
set(gca,'ydir','reverse');


function [enum] = enum_nodos(nelx,nely)
%
%        |   1   |   2   |   3   |   4   |
% enum = |  elx  |  ely  | x.dof | y.dof |
%        |       |       |       |       |
%        |       |       |       |       |
% 
nodo=zeros(((nelx+1)*(nely+1)),4);
c=0;
for ely = 0:nely
    for elx = 0:nelx
        c=c+1;
        nodo(c,1)=elx;
        nodo(c,2)=ely;
    end
end
c=0;
for ely=1:1:nely
    for elx=1:1:nelx
        
        n1 = (nely+1)*(elx-1)+ely; 
        n2 = (nely+1)* elx   +ely;
        c=c+1;
        nodo(c,3)=2*n1-1;
        nodo(c,4)=2*n1  ;
        
        if elx==nelx
            c=c+1;
            nodo(c,3)=2*n2-1;
            nodo(c,4)=2*n2  ;
        end
    end
    
    if ely==nely
        for elx=1:1:nelx
            n1 = (nely+1)*(elx-1)+ely; 
            n2 = (nely+1)* elx   +ely;
            c=c+1;
            nodo(c,3)=2*n1+1;
            nodo(c,4)=2*n1+2;
                
            if elx==nelx
                c=c+1;
                nodo(c,3)=2*n2+1; 
                nodo(c,4)=2*n2+2;
            end
        end
    end 
end
enum=nodo;