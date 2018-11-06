function [Svm] = computestress(U,edofMat,E0,nu,penal,xPhys)
%----------------------------------------------------------
L = [ 1 0 0 0; 0 0 0 1; 0 1 1 0];
% gradient of shape function (isoparametric)
dNx = 1/4*[-1  1 1 -1];
dNy = 1/4*[-1 -1 1  1];
% jacobian
J = [dNx;dNy]*[-1/2 1/2 1/2 -1/2;-1/2 -1/2 1/2 1/2]';
% inverse of the jacobian    
detJ = J(1,1)*J(2,2) - J(1,2)*J(2,1);
invJ = 1/detJ*[J(2,2) -J(1,2); -J(2,1) J(1,1)];
G = [invJ zeros(2); zeros(2) invJ];
% gradient of shape function for isoparametric elements
dN(1,1:2:8) = dNx; dN(2,1:2:8) = dNy;
dN(3,2:2:8) = dNx; dN(4,2:2:8) = dNy; 
% B matrix
B = L*G*dN;
% Von Mises Matrix
V   = [1 -1/2 0; -1/2 1 0; 0 0 3];
D   = [E0/(1-nu*nu),     E0*nu/(1-nu*nu),    0;
       E0*nu/(1-nu*nu),  E0/(1-nu*nu),       0;
       0,                0,               E0/2/(1+nu)];
T   = D*B;
M   = T'*V*T;
%-------------------------------------------------------
% Equivalent stresses
SEQ    = (sum((U(edofMat)*M).*U(edofMat),2)).^0.5;
Svm    = xPhys;
Svm(:) = SEQ(:).*xPhys(:).^penal;
%-------------------------------------------------------
subplot(2,1,1)
colormap(gca,'gray'); imagesc(1-xPhys); colorbar; caxis([0 1]); axis equal; axis off; drawnow;
subplot(2,1,2)
colormap(gca,'jet'); imagesc(Svm); colorbar; axis equal; axis off; drawnow;