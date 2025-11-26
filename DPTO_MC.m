function DPTO_MC(lelx, lely, vol_con, loadcase)
%% PARAMETERS DEFINITION
scale = 1; h = 1;
lelx = 160*scale; lely = 100*scale; loadcase = 1;
BDY = [0, 0; lelx, lely] ;
vol_con = 0.5;
xmin = 1.e-9; dbeta = 1;
E0 = 10; Emin = 1e-8; nu = 0.3;
rmin =  max(1.2*h, scale * 10 * vol_con); % FILTER RADIUS
maxedge = 20*h; minedge = 2*h;
d1 = 0.5*h; d2 = 1.0*h;
%% INITIAL MESHING
[xn, yn] = meshgrid(0: h :lelx, 0: h: lely) ;
[p, t, ~, ~, Ve, pmid] = GenerateMesh(xn, yn, h, BDY,[]) ;
% INITIAL DESIGN VARIABLES
xupp = ones(length(t), 1); xlow = ones(length(t), 1)*xmin;
x = max(xlow,min(xupp,ones(length(t),1)*vol_con));
clf ; colormap summer ; patch( 'Faces', t,'Vertices', p, 'FaceVertexCData', x, 'FaceColor', 'flat') ; colorbar ;
[ H, Hs] = FilterIndex(pmid, rmin, Ve);  % PREPARE FILTER
%% MAIN LOOP
loop = 0 ; loopbeta = 0; beta = 1e-6; tol = 0;
while tol == 0
    xold = x; loop = loop + 1 ; loopbeta = loopbeta + 1;
    % FEA AND SENSITIVITY ANALYSIS
    [dCe, J0] = FEA(t, p, BDY, x, E0, Emin, nu, loadcase) ;
    dVe = Ve./sum(Ve); % SENSITIVITY OF VOLUME
    vol = sum( Ve.* x)/( lelx * lely) ;
    % GRAYNESS DEGREE
    WetE = zeros(size(x));
    for i = 1:length(x)
        WetE(i) = Ve(i).* min(abs(x(i) - 0), abs(x(i) - 1));
    end
    delta = sum(WetE)/sum(Ve);
    % OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES
    [x] = OcUpdate(xold, dCe, dVe, Ve, vol_con, lelx, lely,H,Hs,beta) ;
    % RESULT DISPLAY
    clf; colormap summer ;
    patch( 'Faces', t, 'Vertices', p, 'FaceVertexCData', x, 'FaceColor', 'flat') ;
    colorbar ; hold on ; axis off equal tight ; drawnow ;
    change = max(abs(xold(:)-x(:))) ;
    fprintf( 'It.:%5i Obj.:%8.4f Vol.:%7.3f ch.:%7.3f 0/1:%7.3f\n', loop, J0, vol, change, delta) ;
    % INCREASE PROJECTION PARAMETER
    if ((loopbeta >= 5 && change < 0.02) || loopbeta >= 50)
        if delta < 0.02 && change < 0.02
            tol = 1;
        end
        if delta >= 0.02
            beta = beta + dbeta; loopbeta = 0;
            fprintf('Parameter beta increased to %g.\n',beta);
            if delta < 0.1
                xBF = griddata(pmid(:,1), pmid(:,2), x , xn, yn, 'linear') ;
                xBD = griddata(pmid(:,1), pmid(:,2), x , xn, yn, 'nearest') ;
                xBF(isnan(xBF)) = xBD(isnan(xBF));
                [c] = ContourPoints(contour(xn, yn, xBF, [0.5 0.5]), d1, d2) ;
                [p,t,t1,t2,Ve,pmid] = GenerateMesh(xn, yn, h, BDY, c, xBF, maxedge, minedge) ;    % REMESHING THE DESIGN DOMAIN
                patch('Faces', t1, 'Vertices', p, 'EdgeColor', 'k', 'FaceColor',[0 127 102]/255) ;  axis off equal tight
                patch('Faces', t2, 'Vertices', p, 'EdgeColor', 'k', 'FaceColor',[255 255 102]/255) ; set(gca,'YDir','reverse') ;
                x = ones(length(t), 1); x(1:length(t1)) = xmin ;
                [ H, Hs] = FilterIndex(pmid, rmin, Ve); % PREPARE FILTER FOR NEW MESH
            end
        end
    end
end
%% FINAL RESULT PLOTTTING
xBF = griddata(pmid(:,1), pmid(:,2), x , xn, yn, 'linear') ;
xBD = griddata(pmid(:,1), pmid(:,2), x , xn, yn, 'nearest') ;
xBF(isnan(xBF)) = xBD(isnan(xBF));
[c] = ContourPoints(contour(xn, yn, xBF, [0.5 0.5]), d1, d2) ;
[p,~,t1,t2] = GenerateMesh(xn, yn, h, BDY, c, xBF, maxedge, minedge) ;
patch('Faces', t1, 'Vertices', p, 'EdgeColor', 'none', 'FaceColor',[0 127 102]/255) ;  axis off equal tight
patch('Faces', t2, 'Vertices', p, 'EdgeColor', 'none', 'FaceColor',[255 255 102]/255) ; set(gca,'YDir','reverse') ;
end

%% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES
function [xnew] = OcUpdate(xold, dCe, dVe, Ve, vol_con, lelx, lely, H, Hs, beta)
l1 = 0 ; l2 = 1e10 ;  move = 0.2 ;
while (l2-l1)/(l1+l2) > 1.0e-6
    lmid = 0.5*(l2+l1) ;
    x = max( 0, max( xold-move, min( 1, min( xold+move, xold.*sqrt(max(0,-lmid.*dCe./dVe)))))) ;
    xfilt = max(0,min(1,H*x(:)./Hs));
    xnew = (tanh(beta*0.5)+tanh(beta*(xfilt-0.5)))/2/tanh(beta*0.5);
    xnew = max( 0, max(xold-move, min( 1, min(xold+move, xnew)))) ;
    if sum(xnew(:).*Ve(:)) < vol_con*lelx*lely, l1 = lmid; else l2 = lmid ; end
end
end

%% FIND POINTS ON THE CONTOUR
function [c] = ContourPoints(c,d1,d2)
num = 1; col = 1;
while col < size(c,2)
    idx = col+1:col+c(2,col);
    s(num).ck = c(:,idx);
    s(num).sed = norm(s(num).ck(:,1)-s(num).ck(:,end));
    s(num).isopen = s(num).sed>1e-12;
    num = num+1; col = col+c(2,col)+1;
end
c = [];
for k = 1:num-1
    ck = s(k).ck;
    if length(ck)>4
        ndist = vecnorm(diff(ck,1,2),2,1);
        ck(:,find(ndist < d1) +1) = 0.5*(ck(:,ndist < d1)+ck(:,find(ndist < d1) +1));
        ck(:,ndist < d1) = [];
        if s(k).sed < d1
            ck(:, end) = []; % REMOVE TOO-CLOSE NODES
        end
        if ~s(k).isopen
            ck = [ck ck(:,1)];
        end
        ndist = vecnorm(diff(ck,1,2),2,1);
        ct = ck; ck = ct(:,1);
        for i = 1:length(ndist)
            if  ndist(i) > d2
                ck = [ck 0.5*(ct(:,i)+ct(:,(i+1))) ct(:,i+1)]; % INCLUDE ADDITIONAL NODES
            else
                ck = [ck ct(:,i+1)];
            end
        end
        c = [c; ck(:,1:end-(~s(k).isopen))'];
    end
end
end

%% BODY FITTED MESH GENERATOR
function [p,t,t1,t2,Ve,pmid] = GenerateMesh(xn, yn, h, BDY, c, dN,maxedge,minedge)
if isempty(c) == 1
    [x,y] = meshgrid(BDY(1,1):h: BDY(2,1), BDY(1,2):h: BDY(2,2)) ;
    x(2:2:end,:)=x(2:2:end,:)+0.5*h; pi=[x(:),y(:)];
    pi(:,1) = min(BDY(2,1),max(BDY(1,1),pi(:,1))); pi(:,2) = min(BDY(2,2),max(BDY(1,2),pi(:,2)));
    p = [ BDY(2,1) BDY(1,2); BDY(1,1) BDY(2,2); BDY(1,1) BDY(1,2); BDY(2,1) BDY(2,2); BDY(2,1),(BDY(1,2)+BDY(2,2))/2; pi];
    p = unique( p, 'rows', 'stable');
    t = delaunayn(p);
    t1 = [] ; t2 = [] ;
else
    [x,y] = meshgrid(BDY(1,1): h: BDY(2,1), BDY(1,2): sqrt(3)/2*h : BDY(2,2)) ;
    x(2:2:end,:)=x(2:2:end,:)+0.5*h; pi=[x(:),y(:)];
    d = zeros(size(pi,1),1);
    for i = 1:size(pi,1)
        d(i) = sqrt(min((pi(i,1)-c(:,1)).^2+(pi(i,2)-c(:,2)).^2));
    end
    r0 = 1./min(max(minedge,d),maxedge).^2;
    pfix=[c; BDY(2,1) BDY(1,2); BDY(1,1) BDY(2,2); BDY(1,1) BDY(1,2); BDY(2,1) BDY(2,2); BDY(2,1),(BDY(1,2)+BDY(2,2))/2];
    pfix = unique(pfix,'rows','stable');
    load('randommatrixmc.mat','ranmatrix');   % LOAD SAVED RANDOM MATRIX
    % ranmatrix = rand(size(pi,1),1);       % OR GENERATE A NEW RANDOM MATRIX
    p=[pfix; pi(ranmatrix<r0./max(r0),:)];  % DELETE NODES USING REJECTION METHOD
    p1 = 0;
    % NOVE-MOVING LOOPS
    for i = 1:500
        if max(sum((p-p1).^2,2))>0.1*h
            p = unique(p,'rows','stable');
            t = delaunayn(p);  % DELAUNAY TRIANGULATION
            edges = unique(sort([t(:,[1,2]);t(:,[1,3]);t(:,[2,3])],2),'rows') ;
            p1 = p;
            midpoint = (p(edges(:,1),:)+p(edges(:,2),:))/2;
            d = zeros(size(midpoint,1),1);
            for j = 1:size(midpoint,1)
                dist = (midpoint(j, 1) - c(:, 1)).^2 + (midpoint(j, 2) - c(:, 2)).^2;
                d(j) = sqrt(min(dist));
            end
            L1 = min(max(minedge,d),maxedge);
        end
        L = sqrt(sum((p(edges(:,1),:)-p(edges(:,2),:)).^2,2));
        L0 = 1.2*L1*sqrt(sum(L.^2)/sum(L1.^2));
        Fb = max(L0-L,0)./L *[1,1].*(p(edges(:,1),:)-p(edges(:,2),:));
        Fp = full(sparse(edges(:,[1,1,2,2]),ones(size(d))*[1,2,1,2],[Fb,-Fb],size(p,1),2));
        Fp(1:size(pfix,1),:) = 0 ;
        p = p + 0.2*Fp ;  % MOVE NODE ACCORDING TO VIRTUAL FORCE
        p(:,1) = min(BDY(2,1),max(BDY(1,1),p(:,1))) ;
        p(:,2) = min(BDY(2,2),max(BDY(1,2),p(:,2))) ;
        if max(sum(0.2*Fp,2))<0.05*h, break; end
    end
    pmid = (p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;
    dE = interp2(xn,yn,dN,pmid(:,1),pmid(:,2),'cubic');  % INTERPOLATE DENSITY INTO ELEMENT CENTROIDS
    t1=[t(dE<0.5,:)]; t2=t(dE>=0.5,:); t=[t1;t2];
end
pmid = (p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;
Ve = zeros(length(t),1) ;
for j = 1:length(t)
    Ve(j) = 0.5.*det([ones(3,1) p(t(j,:),:)]);  % ELEMENT VOLUME
end
end

%% PREPARE FILTER
function [H, Hs] = FilterIndex(pmid, rmin,Ve)
maxNum =  ceil ( 100 * rmin);
iH = ones(length(pmid)*maxNum, 1) ;
jH = ones(length(pmid)*maxNum, 1) ;
sH = zeros(length(pmid)*maxNum, 1) ;
Count = 1;
for i = 1 : length(pmid)
    pi = pmid(i,:);
    pni = find(abs(pmid(:,1)-pi(1))<=rmin & abs(pmid(:,2)-pi(2))<=rmin);
    Weight = max((rmin - pdist2(pi,pmid(pni,:))),0) .* (Ve(pni))';
    Countnew = Count + length(Weight);
    iH(Count:Countnew-1) = i;
    jH(Count:Countnew-1) = pni;
    sH(Count:Countnew-1) = Weight;
    Count = Countnew;
end
H = sparse( iH, jH, sH) ;
Hs = sum( H, 2) ;
end

%% ELEMENT STIFFNESS MATRIX
function [Ke] = ElementMatrixKe(X,Y,E0,nu)
D = E0/(1-nu^2)*[1  nu  0;  nu  1  0;  0  0  (1-nu)/2] ;
J = [X(1)-X(3) Y(1)-Y(3) ; X(2)-X(3) Y(2)-Y(3)] ;
Be = 1/det(J)*[J(2,2) 0 -J(1,2) 0 -J(2,2)+J(1,2) 0; 0 -J(2,1) 0 J(1,1) 0 J(2,1)-J(1,1); -J(2,1) J(2,2) J(1,1) -J(1,2) J(2,1)-J(1,1) -J(2,2)+J(1,2)];
Ke = 1/2*det(J)*Be'*D*Be ;
end

%% FINITE ELEMENT ANALYSIS
function [dCe, J] = FEA(t, p, BDY, x, E0, Emin, nu, loadcase)
NT = length(t); KK = zeros( 6, 6*NT) ;
for pi = 1:NT
    KK(:,6*pi-5:6*pi) = (Emin+x(pi)*(E0-Emin)) * ElementMatrixKe(p(t(pi,:),1),p(t(pi,:),2),1,nu) ;
end
elemDof = zeros(NT,6); elemDof(:,[1 3 5]) = 2*t-1; elemDof(:,[2 4 6]) = 2*t;
iK = reshape( kron( elemDof, ones(6,1))', 36*NT,1) ;
jK = reshape( kron( elemDof, ones(1,6))', 36*NT,1) ;
sK = reshape( KK, 36*NT, 1) ;
NK = sparse( iK, jK, sK, 2*length(p), 2*length(p)) ;
NK = (NK+NK')/2 ;
%% LOAD CASES
if loadcase == 1  % CANTILEVER BEAM (160 X 100)
    fixedNodes = find(p(:,1)==BDY(1,1)) ;
    forceNodes = find(p(:,1)==BDY(2,1) & p(:,2)==(BDY(1,2)+BDY(2,2))/2) ;
    fixedDof = [2*fixedNodes-1; 2*fixedNodes] ;
    AllDof = 1:2*length(p) ;
    freeDofs = setdiff( AllDof, fixedDof) ;
    F = sparse( 2*forceNodes, 1, -1, 2*length(p), 1) ;
elseif loadcase == 2  % MBB BEAM (240 X 80)
    fixedNodes1 = find(p( :, 1) == BDY( 1, 1)) ;
    fixedDof1 = 2*fixedNodes1 - 1 ;
    fixedNodes2 = find(p(:,1) == BDY(2,1) & p(:,2) == BDY(2,2)) ;
    fixedDof2 = 2*fixedNodes2 ;
    fixedDof = [fixedDof1;fixedDof2] ;
    forceNodes = find(p(:,1)==BDY(1,1) & p(:,2)==BDY(2,2)) ;
    AllDof = 1:2*length(p) ;
    freeDofs = setdiff( AllDof, fixedDof) ;
    F = sparse( 2*forceNodes, 1, -0.5, 2*length(p), 1) ;
end
U = zeros( 2*length(p), 1) ;
U( freeDofs, :) = NK(freeDofs,freeDofs) \ F(freeDofs,1) ;
for pi = 1 : NT
    Ce(pi) = sum((U(elemDof(pi,:))'*KK(:,6*pi-5:6*pi)).*U( elemDof(pi,:))', 2) ;
    dCe(pi) = -Ce(pi)./x(pi); % SENSITIVITY OF COMPLIANCE
end
J = 0.5.*sum(Ce); dCe = dCe';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by ZC. Zhuang and Y. M. Xie                %
% Centre for Innovative Structures and Materials, RMIT University         %
% Please send your comments to: zicheng.zhuang@rmit.edu.au                %
%                                                                         %
% The program is proposed for educational purposes and is introduced in   %
% the paper - Improving Stress Distribution in Topology                   %
% Optimization Using Density Projection, Multiple Constraints,            %
% and Body-Fitted Mesh, CMAME, 2025                                       %
%                                                                         %
% Disclaimer:                                                             %
% The authors reserve all rights but do not guarantee that the code is    %
% free from errors. Furthermore, we shall not be liable in any event.     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%