function DPTO_STR(lelx, lely, str_con)
%% PARAMETERS DEFINITION
scale = 1; h = 1;
lelx = 150*scale; lely = 150*scale;  % LENGTH OF DESIGN DOMAIN
lpd = 60*scale; lload = 6*scale;    % LENGTH OF PASSIVE DOMAIN AND LOAD REGION
BDY = [0, 0; lelx, lely] ;
vol_ini = 0.3; str_con = 0.4;   % INITIAL VOLUME AND STRESS CONSTRAINT
xmin = 1e-9; dbeta = 1;  move = 0.02;
E0 = 1; Emin = 1e-9; nu = 0.3; eps = 0.03;
rmin =  2.5 * scale; pnorm = 16;  % FILTER RADIUS AND P-NORM EXPONENT
maxedge = 40*h; minedge = 2*h;
d1 = 0.5*h; d2 = 1.0*h;
%% INITIAL MESHING
[xe, ye] = meshgrid(0.5: h :lelx-0.5, 0.5: h: lely-0.5); % REGULAR ELEMENT CENTROID IN DESIGN DOMAIN
[xn, yn] = meshgrid(0: h :lelx, 0: h: lely); % REGULAR NODES IN DESIGN DOMAIN
[p, t, pv, ps, pd, td, ~, ~, Ve, pmid] = GenerateMesh(xn, yn, h, BDY, lpd, lload, []);
% LOWER AND UPPER BOUNDS OF DESIGN VARIABELS
xupp = ones(length(t), 1); xlow = ones(length(t), 1)*xmin;
xlow(ps) = 1; xupp(pv) = xmin;
% INITIAL DESIGN VARIABLES
x = max(xlow,min(xupp,ones(length(t),1)*vol_ini));
xphys = x; xfilt = x;
[ H, Hs] = FilterIndex(pmid, rmin, Ve ,xupp, xmin); % PREPARE FILTER
%% MAIN LOOP
loop = 0 ; loopbeta = 0; loopcon = 0; beta = 1e-6; tol = 0; remesh = 0; loopremesh = 0; U = zeros(2*length(p),3);
%% INITIALIZE MMA OPTIMIZER
m       = 2;                        % NUMBER OF CONSTRAINTS
n       = length(td);               % NUMBER OF ELEMENTS
xdold   = x(pd);                    % OLD DESIGN VARIABLE
a0      = 1;
ai      = zeros(m,1);
ci      = 1000*ones(m,1);
di      = zeros(m,1);
low     = zeros(n,1);
upp     = ones(n,1);
while tol == 0
    xold = x; oldxp = xphys;       % OLD DENSITY VALUES
    loop = loop + 1 ; loopbeta = loopbeta + 1;  loopcon = loopcon + 1;  beta1 = 1.2*beta;
    %% FEA AND SENSITIVITY ANALYSIS
    [U0, Compliance, str1, vms, S, K, elemDof] = FEA(t, p, BDY, xphys, E0, Emin, nu, lpd, lload) ;
    U(:,1) = U0(:);
    for e = 1:length(t)
        str_star(e,1) = xphys(e,1)*vms(e,1)./((xphys(e,1)+(1-xphys(e,1))*eps));   % VON MISES STRESS
        star_x(e,1) = eps*vms(e,1)./((xphys(e,1)+(1-xphys(e,1))*eps).^2);   % d VON MISES STRESS/d X
        star_vms(e,1) = xphys(e,1)./((xphys(e,1)+(1-xphys(e,1))*eps));
        star_s1(e,:) = star_vms(e,1)/(vms(e,1))*[str1(e,1)-0.5*str1(e,2) str1(e,2)-0.5*str1(e,1) 3*str1(e,3)]';
        Se{e} = S(:,6*e-5:6*e);
        Ke{e} = K(:,6*e-5:6*e);
    end
    vms_max(loop) = max(str_star);   % MAXIMAL VON MISES STRESS
    vms_pn(loop) = (sum(str_star.^pnorm))^(1/pnorm);   % P-NORM STRESS
    pn_star = (sum(str_star.^pnorm))^(1/pnorm-1)*(str_star.^(pnorm-1));
    pn_x = 0;
    vms_aver(loop) = sum(str_star.*Ve)./sum(xphys.*Ve); % OVERALL VON MISES STRESS
    aver_star = Ve./sum(xphys.*Ve);
    aver_x = -Ve.*sum(str_star.*Ve)./(sum(xphys.*Ve))^2;
    %% STRESS VALUE INTERPOLATION
    if remesh == 0
        strBF = griddata(pmid(:,1), pmid(:,2), str_star , xe, ye, 'cubic') ;
    else
        strBF = griddata(pmid(:,1), pmid(:,2), str_star , xe, ye, 'linear') ;
    end
    strBD = griddata(pmid(:,1), pmid(:,2), str_star , xe, ye, 'nearest') ;   % STRESS INTERPOLATION FOR COMPARISON
    strBF(isnan(strBF)) = strBD(isnan(strBF));
    strBF_aver = sum(sum(strBF))/sum(xphys.*Ve); % OVERALL VON MISES STRESS
    vms_max_rec = max(max(strBF));  % MAXIMAL VON MISES STRESS IN RECTANGLES
    figure(2); clf; imagesc(strBF);colormap(gca,'jet'); set(gca,'YDir','normal');  % PLOT STRESS DISTRIBUTION
    title(['max stress = ',num2str(max(strBF(:))),';vmsav = ',num2str(strBF_aver)]);
    colorbar;  axis equal off; pause(1e-6);
    %% SOLVING ADJOINT STATES
    F1(1:2*length(p),1) = 0; F2(1:2*length(p),1) = 0; % GLOBAL AND ELEMENT DUMMY FORCE
    nodF1(1:length(t),1:6) = 0;   nodF2(1:length(t),1:6) = 0;
    for e = 1:length(t)
        nodF1(e,:) = pn_star(e,1)*star_s1(e,:)*Se{e};
        F1(elemDof(e,:),1) = F1(elemDof(e,:),1)+[nodF1(e,1);nodF1(e,2);nodF1(e,3);nodF1(e,4);nodF1(e,5);nodF1(e,6)];
        nodF2(e,:) = aver_star(e,1)*star_s1(e,:)*Se{e};
        F2(elemDof(e,:),1) = F2(elemDof(e,:),1)+[nodF2(e,1);nodF2(e,2);nodF2(e,3);nodF2(e,4);nodF2(e,5);nodF2(e,6)];
    end
    [U1] = FEA_AJ(t, p, BDY, xphys, F1, E0, nu, lpd);
    U(:,2) = U1(:);
    [U2] = FEA_AJ(t, p, BDY, xphys, F2, E0, nu, lpd);
    U(:,3) = U2(:);
    %% DERAVATIVES FOR VOLUME And TWO STRESS CONSTRAINTS
    dpn(1:length(t),1) = 0; dav(1:length(t),1) = 0;
    for e = 1:length(t)
        Ue = U(elemDof(e,:),1);
        Ue1 = U(elemDof(e,:),2);
        Ue2 = U(elemDof(e,:),3);
        dpn(e,1) = pn_x + pn_star(e,1)*star_x(e,1)-Ue1'*Ke{e}*Ue;
        dav(e,1) = aver_x(e,1) + aver_star(e,1) *star_x(e,1)-Ue2'*Ke{e}*Ue;
    end
    dv = Ve./sum(xupp.*Ve);
    %% OBJECTIVE FUNCTION AND GRAYNESS DEGREE
    vol(loop) = sum(Ve.* xphys)/ sum(xupp.*Ve) ;
    WetE = zeros(size(xphys));
    for i = 1:length(xphys)
        WetE(i) = Ve(i).* min(abs(xphys(i) - 0), abs(xphys(i) - 1));
    end
    delta = sum(WetE)/sum(Ve);  % CALCULATE GREYNESS DEGREE
    %% SENSITIVITY TRANSFORMATION
    dv(:) = dv(:).*(1-tanh(beta1*(xfilt(:)-0.5)).^2)*beta1/2/tanh(beta1*0.5);
    dpn(:) = dpn(:).*(1-tanh(beta1*(xfilt(:)-0.5)).^2)*beta1/2/tanh(beta1*0.5);
    dav(:) = dav(:).*(1-tanh(beta1*(xfilt(:)-0.5)).^2)*beta1/2/tanh(beta1*0.5);
    dv = H*dv(:)./Hs;
    dpn = H*dpn(:)./Hs;
    dav = H*dav(:)./Hs;
    %% MMA UPDATING OF DESIGN VARIABLES
    xval = x(pd);
    dval = dv(pd); dpnal = dpn(pd); daval = dav(pd);
    f0val=vol(loop)./vol(1+remesh*loopremesh);  % DEFINE OBJECTIVE FUNCTION VALUE
    df0dx=dval./vol(1+remesh*loopremesh);
    fval = zeros(m,1);
    if (mod(loopcon,10) == 1) || (loopbeta ==1)
        g1_con = vms_pn(loop).*(str_con /vms_max_rec);   % DEFINE P-NORM STRESS CONSTRAINT VALUE
        if remesh == 0
            g2_con = vms_max_rec;  % DEFINE OVERALL STRESS CONSTRAINT VALUE
        else
            g2_con = 0.42;  % LARGER STRESS VALUE AFTER REMESH
        end
    end
    fval(1,1) = vms_pn(loop)./g1_con-1;
    dfdx(1,1:n) = dpnal(:)'./g1_con;
    fval(2,1) = vms_aver(loop)./g2_con-1;
    dfdx(2,1:n) = daval(:)'./g2_con;
    xxmax = ones(n,1); xxmin = ones(n,1)*xmin;
    [xmma,~,~,~,~,~,~,~,~,low,upp] = mmasub(m,n,loop-remesh*loopremesh,xval,xxmin,xxmax,xdold,xdold,f0val,df0dx,fval,dfdx,low,upp,a0,ai,ci,di);
    x(pd) = xmma;
    %% FILTERING AND PROJECTION
    x=max(xlow,min(xupp,x));
    xfilt = max(xlow,min(xupp,H*x(:)./Hs));  % DENSITY FILTERING
    xphys = (tanh(beta1*0.5)+tanh(beta1*(xfilt-0.5)))/2/tanh(beta1*0.5);  % DENSITY PROJECTION
    xphys = max(xlow,min(xupp,max(oldxp-move,min(oldxp+move,xphys))));
    l1 = 0; l2 = 1; sum_x = sum(x(:));
    while (l2-l1) > 1.0e-10
        th = (l1+l2)/2.0;
        x = max(xlow,min(xupp,max(xold-move,min(xold+move,(tanh(beta*th)...
            +tanh(beta*(xfilt-th)))/(tanh(beta*th)+tanh(beta*(1.-th)))))));
        if sum(x(:)) > sum_x; l1 = th; else; l2 = th; end   %  KEEP DESIGN VARIABLE VOLUME CONSTANT
    end
    xdold=x(pd);
    figure(1); clf; colormap summer ;
    patch('Faces', t, 'Vertices', p, 'EdgeColor', 'none', 'FaceVertexCData', xphys, 'FaceColor', 'flat') ; % PLOT OPTIMIZED STRUCTURE
    colorbar;  axis off equal tight ; drawnow ;
    fprintf('It.:%5i Obj.:%8.4f Vol.:%7.3f vmmax.:%7.5f vmpn.:%7.5f vmav.:%7.5f 0/1:%7.5f beta:%5.2f\n', ...
        loop, Compliance, vol(loop), vms_max_rec, vms_pn(loop), vms_aver(loop), delta, beta) ;
    %% CONVERGENCE AND TERMINATION CHECK
    if loopbeta > 10
        obj_ch = abs(mean(vol(loop-9:loop-5))-mean(vol(loop-4:loop)))./mean(vol(loop-9:loop));
        cons_ch = abs(mean(vms_max(loop-9:loop-5))-mean(vms_max(loop-4:loop)))./mean(vms_max(loop-9:loop));
        vms_max_ch = abs(vms_max_rec-str_con);
        if remesh == 0
            %% INCREASING PROJECTION PARAMETER WHEN CONVERGED
            if delta >= 0.02 && cons_ch < 0.005 && obj_ch < 0.001
                beta = beta + dbeta;
                loopbeta = 0;
                fprintf('Parameter beta increased to %g.\n',beta);
            end
            if delta >= 0.02 && loopbeta >= 50
                beta = beta + 0.5 * dbeta;
                loopbeta = 0;
                fprintf('Parameter beta increased to %g.\n',beta);
            end
            %% REMESHING THE DESIGN DOMAIN
            if delta < 0.02 && cons_ch < 0.005 && obj_ch < 0.001 && vms_max_ch < 0.01   % WHETHER OR NOT TO REMESH
                remesh = 1;
                xBF = griddata(pmid(:,1), pmid(:,2), x , xe, ye, 'linear') ;
                xBD = griddata(pmid(:,1), pmid(:,2), x , xe, ye, 'nearest') ;
                xBF(isnan(xBF)) = xBD(isnan(xBF));
                [xBF] = Expand(xBF);
                [c] = ContourPoints(contour(xn, yn, xBF, [0.5 0.5]), d1, d2) ;
                [p,t, pv, ps, pd, td, t1,t2,Ve,pmid,dE] = GenerateMesh( xn, yn, h, BDY, lpd, lload, c, xBF,maxedge,minedge);
                % REDINFING LOWER AND UPPER BOUNDS OF DESIGN VARIABELS
                xupp = ones(length(t), 1); xlow = ones(length(t), 1)*xmin;
                xlow(ps) = 1;
                xupp(pv) = xmin;
                x = ones(length(t), 1);
                x(dE<0.5) = xmin;
                x = max(xlow,min(xupp,x));
                xphys = x; xfilt = x;
                [ H, Hs] = FilterIndex(pmid, rmin, Ve ,xupp,xmin); % PREPARE FILTER FOR NEW MESH
                U = zeros(2*length(p),3);
                str_star = []; dpn = []; dav = []; star_x=[]; star_vms = []; star_s1 = []; Se = []; Ke = [];
                F1 = []; F2 = []; nodF1=[]; nodF2=[]; dfdx = []; % RESET ALL STRESS FIELDS
                loopbeta = 0; loopcon = 0; loopremesh = loop;
                beta = beta - 0.5 * dbeta;
                %% REINITIALIZE MMA OPTIMIZER
                n       = length(td);
                xdold   = x(pd);
                low     = zeros(n,1);
                upp     = ones(n,1);
            end
        else
            if delta < 0.02 && cons_ch < 0.005 && obj_ch < 0.001 && vms_max_ch < 0.005  % TERMINATION CHECK (STRICTER)
                tol = 1;
            end
        end
    end
end
%% FINAL RESULT PLOTTTING
xBF = griddata(pmid(:,1), pmid(:,2), x , xe, ye, 'linear') ;
xBD = griddata(pmid(:,1), pmid(:,2), x , xe, ye, 'nearest') ;
xBF(isnan(xBF)) = xBD(isnan(xBF));
[xBF] = Expand(xBF);
[c] = ContourPoints(contour(xn, yn, xBF, [0.5 0.5]), d1, d2) ;
[p, ~, ~, ~, ~, ~, t1,t2] = GenerateMesh( xn, yn, h, BDY, lpd, lload, c, xBF,maxedge,minedge);
figure(3); clf; patch('Faces', t1, 'Vertices', p, 'EdgeColor', 'none', 'FaceColor',[0 127 102]/255) ;  axis off equal tight
patch('Faces', t2, 'Vertices', p, 'EdgeColor', 'none', 'FaceColor',[255 255 102]/255) ;
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
function [p,t,pv,ps,pd,td,t1,t2,Ve,pmid,dE] = GenerateMesh(xn,yn,h,BDY,lpd,lload,c,dN,maxedge,minedge)
if isempty(c) == 1
    [x,y] = meshgrid(BDY(1,1):h: BDY(2,1), BDY(1,2):h: BDY(2,2));
    pi=[x(:),y(:)];
    pi(:,1) = min(BDY(2,1),max(BDY(1,1),pi(:,1))); pi(:,2) = min(BDY(2,2),max(BDY(1,2),pi(:,2)));
    p = [BDY(2,1) BDY(1,2); BDY(1,1) BDY(2,2); BDY(1,1) BDY(1,2); BDY(2,1) BDY(2,2); pi];
    p = unique( p, 'rows', 'stable');
    t = delaunayn(p);
    pmid = (p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;
    pv = pmid(:,1)>lpd & pmid(:,2)>lpd;
    ps = pmid(:,1)>BDY(2,2)-lload & pmid(:,2)>lpd-1 & pmid(:,2)<lpd;
    pd = ~(pv|ps);
    tpv = t((pmid(:,1)>lpd & pmid(:,2)>lpd),:);
    tps = t((pmid(:,1)>BDY(2,2)-lload & pmid(:,2)>lpd-1 & pmid(:,2)<lpd),:);
    td = setdiff(t,union(tpv,tps,'rows'),'rows');
    t1 = [] ; t2 = [] ; dE = [];
else
    [xp1,yp1] = meshgrid(lpd,lpd:BDY(2,2));
    [xp2,yp2] = meshgrid(lpd:BDY(2,1),lpd);
    [xp3,yp3] = meshgrid(BDY(2,2)-lload,lpd-1:lpd);
    [xp4,yp4] = meshgrid(BDY(2,2)-lload:BDY(2,2),lpd-1);
    fixed_passive = [xp1(:),yp1(:); xp2(:),yp2(:); xp3(:),yp3(:); xp4(:),yp4(:)];
    fixed_passive = unique(fixed_passive, 'rows', 'stable');  % FIX NODES ON THE BOUNDARIES
    [x,y] = meshgrid(BDY(1,1): h: BDY(2,1), BDY(1,2): sqrt(3)/2*h : BDY(2,2)) ;
    c1 = c (c(:,1)>=59 & c(:,1)<=61 & c(:,2)>=59.8,:);
    c2 = c (c(:,2)>=59 & c(:,2)<=61 & c(:,1)>=140,:);
    c = setdiff(c,c1,'rows');
    c = setdiff(c,c2,'rows');   % REMOVE CONTOUR NODES TOO CLOSE TO THE BOUNDARIES
    c = [fixed_passive; c];
    c = unique(c,'rows','stable');
    x(2:2:end,:)=x(2:2:end,:)+0.5*h; pi=[x(:),y(:)];
    d = zeros(size(pi,1),1);
    for i = 1:size(pi,1)
        d(i) = sqrt(min((pi(i,1)-c(:,1)).^2+(pi(i,2)-c(:,2)).^2));
    end
    r0 = 1./min(max(minedge,d),maxedge).^2;
    pfix=[fixed_passive; c; BDY(2,1) BDY(1,2); BDY(1,1) BDY(2,2); BDY(1,1) BDY(1,2); BDY(2,1) BDY(2,2)];
    pfix = unique(pfix,'rows','stable');
    load('randommatrix.mat','ranmatrix');   % LOAD SAVED RANDOM MATRIX
    % ranmatrix = rand(size(pi,1),1);       % OR GENERATE A NEW RANDOM MATRIX
    p=[pfix; pi(ranmatrix<r0./max(r0),:)];  % DELETE NODES USING REJECTION METHOD
    p1 = 0;
    % NOVE-MOVING LOOPS
    for i = 1:200
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
        gap = mean(sqrt(sum(0.2*Fp.^2,2)))./h;
        if gap<0.01, break; end
    end
    pmid = (p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;
    pv = pmid(:,1)>lpd & pmid(:,2)>lpd;
    ps = pmid(:,1)>BDY(2,2)-lload & pmid(:,2)>lpd-1 & pmid(:,2)<lpd;
    pd = ~(pv|ps);
    tpv = t((pmid(:,1)>lpd & pmid(:,2)>lpd),:);
    tps = t((pmid(:,1)>BDY(2,2)-lload & pmid(:,2)>lpd-1 & pmid(:,2)<lpd),:);
    td = setdiff(t,union(tpv,tps,'rows'),'rows');
    dE = interp2(xn,yn,dN,pmid(:,1),pmid(:,2),'cubic');  % INTERPOLATE DENSITY INTO ELEMENT CENTROIDS
    t1 = setdiff(union(t(dE<0.5,:),t(pv,:),'rows'),t(ps,:),'rows');  % VOID ELEMENTS
    t2 = setdiff(t,t1,'rows');  % SOLID ELEMENTS
end
Ve = zeros(length(t),1) ;
for j = 1:length(t)
    Ve(j) = 0.5.*det([ones(3,1) p(t(j,:),:)]);  % ELEMENT VOLUME
end
end

%% EXPAND DENSITY FIELD
function [x] = Expand(xs)
xl=[xs(1,1) xs(1,:) xs(1,end); xs(:,1) xs xs(:,end); xs(end,1) xs(end,:) xs(end,end)];
x = 0.25*(xl(1:end-1,1:end-1)+xl(2:end,1:end-1)+xl(1:end-1,2:end)+xl(2:end,2:end));
end

%% PREPARE FILTER
function [H,Hs] = FilterIndex(pmid,rmin,Ve,xupp,xmin)
maxNum =  ceil ( 100 * rmin);
iH = ones(length(pmid)*maxNum, 1) ;
jH = ones(length(pmid)*maxNum, 1) ;
sH = zeros(length(pmid)*maxNum, 1) ;
Count = 1;
for i = 1 : length(pmid)
    pi = pmid(i,:);
    pni = find(abs(pmid(:,1)-pi(1))<=rmin & abs(pmid(:,2)-pi(2))<=rmin);
    Weight = max((rmin - pdist2(pi,pmid(pni,:))),0) .* (Ve(pni))';
    Weight(xupp(pni) == xmin) = 0;
    Countnew = Count + length(Weight);
    iH(Count:Countnew-1) = i;
    jH(Count:Countnew-1) = pni;
    sH(Count:Countnew-1) = Weight;
    Count = Countnew;
end
H = sparse( iH, jH, sH) ;
Hs = sum( H, 2) ;
end

%% ELEMENT MATRIX DERIVATION
function [S,D,Be,Ke] = ElementMatrixKe(X,Y,E0,nu)
D = E0/(1-nu^2)*[1  nu  0;  nu  1  0;  0  0  (1-nu)/2] ;  % ELASTICITY MATRIX
J = [X(1)-X(3) Y(1)-Y(3) ; X(2)-X(3) Y(2)-Y(3)] ;
Be = 1/det(J)*[J(2,2) 0 -J(1,2) 0 -J(2,2)+J(1,2) 0;
    0 -J(2,1) 0 J(1,1) 0 J(2,1)-J(1,1);
    -J(2,1) J(2,2) J(1,1) -J(1,2) J(2,1)-J(1,1) -J(2,2)+J(1,2)];  % STRAIN MATRIX
Ke = 1/2*det(J)*Be'*D*Be;   % STIFFNESS MATRIX
S = D * Be;  % STRESS MATRIX
end

%% FINITE ELEMENT ANALYSIS
function [U,c,str1,vms,S,K,elemDof] = FEA(t,p,BDY,x,E0,Emin,nu,lpd,lload)
NT = length(t); KK = zeros( 6, 6*NT) ;
for i = 1:NT
    [Se,~,~,Ke]  = ElementMatrixKe(p(t(i,:),1),p(t(i,:),2),E0,nu);
    K(:,6*i-5:6*i) = Ke;
    S(:,6*i-5:6*i) = Se;
    KK(:,6*i-5:6*i) = x(i) * Ke ;
end
elemDof = zeros(NT,6); elemDof(:,[1 3 5]) = 2*t-1; elemDof(:,[2 4 6]) = 2*t;
iK = reshape( kron( elemDof, ones(6,1))', 36*NT,1) ;
jK = reshape( kron( elemDof, ones(1,6))', 36*NT,1) ;
sK = reshape( KK, 36*NT, 1) ;
NK = sparse( iK, jK, sK, 2*length(p), 2*length(p)) ;
NK = (NK+NK')/2 ;
%% LOAD CASE (L-BRACKET)
fixedNodes = find(p(:,1)<= lpd & p(:,2)==BDY(2,2)) ;
fixedDof = [2*fixedNodes-1; 2*fixedNodes] ;
AllDof = 1:2*length(p) ;
freeDofs = setdiff( AllDof, fixedDof) ;
forceNodes = find(p(:,1)<BDY(2,1) & p(:,1)>BDY(2,1)-lload & p(:,2)==lpd) ;
forceNodes1 = find((p(:,1)==BDY(2,1) | p(:,1)==BDY(2,1)-lload) & p(:,2)==lpd) ;
F = sparse( 2*length(p),1);
F0 = 1;
F(2*forceNodes) = -F0/lload;
F(2*forceNodes1) = -0.5*F0/lload;
str1 = zeros(NT,3);
vms = zeros(NT,1);
U = zeros( 2*length(p), 1) ;
U( freeDofs, :) = NK(freeDofs,freeDofs) \ F(freeDofs,1) ;
for i = 1 : NT
    Ce(i) = sum((U(elemDof(i,:))'*KK(:,6*i-5:6*i)).*U(elemDof(i,:))', 2) ;
    str1(i,:) = S(:,6*i-5:6*i) * U(elemDof(i,:));
    vms(i,1) = sqrt(str1(i,1)^2+str1(i,2)^2-str1(i,1)*str1(i,2)+3*str1(i,3)^2);
end
c = 0.5.*sum(Ce);
end

%% ADJOINT FEA
function [U] = FEA_AJ(t,p,BDY,x,F,E0,nu,lpd)
NT = length(t); KK = zeros( 6, 6*NT) ;
for i = 1:NT
    [~,~,~,Ke]  = ElementMatrixKe(p(t(i,:),1),p(t(i,:),2),E0,nu);
    KK(:,6*i-5:6*i) = x(i) * Ke ;
end
elemDof = zeros(NT,6); elemDof(:,[1 3 5]) = 2*t-1; elemDof(:,[2 4 6]) = 2*t;
iK = reshape( kron( elemDof, ones(6,1))', 36*NT,1) ;
jK = reshape( kron( elemDof, ones(1,6))', 36*NT,1) ;
sK = reshape( KK, 36*NT, 1) ;
NK = sparse( iK, jK, sK, 2*length(p), 2*length(p)) ;
NK = (NK+NK')/2 ;
fixedNodes = find(p(:,1)< lpd & p(:,2)==BDY(2,2)) ;
fixedDof = [2*fixedNodes-1; 2*fixedNodes] ;
AllDof = 1:2*length(p) ;
freeDofs = setdiff( AllDof, fixedDof) ;
U = zeros( 2*length(p), 1) ;
U( freeDofs, :) = NK(freeDofs,freeDofs) \ F(freeDofs,1) ;
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