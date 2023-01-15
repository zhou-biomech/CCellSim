classdef Cell < handle
    properties
        % nodes = Node.empty;
        p                         % p = [x, y] of nodes
        i, j, k
        x                         % X-axis coordinate of nodes
        y                         % Y-axis coordinate of nodes
        v                         % velocity of nodes
        f                         % nodes' force
        mols = Molecule.empty;    % molecules
        lc                        % Arc length corresponding to p
        dn                        % normal unit vector
        xc = [];                  % centrods for plot
        yc = [];                  % centrods for plot
        
        %% Viscoelasticity of cell membrane
        l, l0                     % equilibrium value of the springs length
        A, A0                     % equilibrium value of the area
        theta0                    % equilibrium value of the angle
        kl                        % spring constant
        ka                        % penalty coefficient of area changes
        kb                        % bending constant
        mu                        %  
        gamma                     % 
        h                         % plot handle
        hc                        % centrod plot handle
        
        %% CSK parameters
        pcsk                      % csk nodes position [x, y]
        ecsk                      % [i, j]   edge node index
        tcsk                      % [i, j, k] triangle node index
        Mcsk                      % transform matrix
        Acsk0                     % equilibrium value of the angle
        Acsk                      % angle
        dAcsk                     % dAcsk = Acsk - Acsk0;
        pAcsk                     % angle variation of points
        hcsk                      % csk plot handle
        rcsk
        Mt2p
        isopen = false;
        res
        res0
        tot_gap
        
        isadh
    end
    
    methods

        function obj = Cell(xc, yc, R, N, totgap)
            [p0, l0, theta0, A0, kl, kb, ka, gamma, mu] = initcell(R, N);
            
            obj.l0 = l0;
            obj.A0 = A0;
            obj.theta0 = theta0;
            obj.kl = kl;
            obj.ka = ka;
            obj.kb = kb;
            obj.mu = mu;
            obj.gamma = gamma;
            obj.p = p0 + [xc, yc];
            obj.v = zeros(size(p0));
            
            obj.isadh = false(length(p0),1);
            
            N = length(p0);
            
            obj.i = [N 1:N-1];       % index of 3 particles on a angle
            obj.j = [1 2:N  ];       % j = i + 1
            obj.k = [2 3:N 1];       % k = j + 1 = i + 2
            
            
            % Generate initial CSK network
            fd = @(p) sqrt(sum(p.^2,2)) - R;
            fh = @(p) R*ones(size(p,1),1);
            dl = 2*pi*R/N;
            
            [pcsk,tcsk] = distmesh2d(fd, fh, dl, R*[-1,-1;1,1], p0);
            
            rcsk = sqrt(sum(pcsk.^2,2)); obj.rcsk = rcsk/max(rcsk);
            
            TR = triangulation(tcsk,pcsk);
            obj.pcsk = pcsk+[xc,yc];
            obj.tcsk = tcsk;
            
            ecsk = edges(TR);
            ecsk(all(ecsk<=N,2), :) = [];
            obj.ecsk = ecsk;

            E = zeros(length(pcsk), length(ecsk)); D = E';
            for k = 1:length(ecsk)
                E(ecsk(k,:),k) = [1, -1]; D(k,ecsk(k,:)) = [-1, 1];
            end
            E(1:N, :) = 0;
            [V, D] = eig(E*D);
            D = D + eye(size(D));
            obj.Mcsk = V*(D==1)/V;
          
            np = length(obj.pcsk);
            nt = length(tcsk);
            Mt2p = sparse(np,nt*3);
            for i = 1:nt
                Mt2p(tcsk(i,1),i) = 1;
                Mt2p(tcsk(i,2),nt+i) = 1;
                Mt2p(tcsk(i,3),nt*2+i) = 1;
            end
            obj.Mt2p = Mt2p;
            
            cellcsk(obj);
            obj.Acsk0 = obj.Acsk;
            obj.pAcsk = zeros(length(obj.pcsk),1);
            
            obj.tot_gap = totgap;
            deg = obj.dAcsk(:);
            obj.res = obj.tot_gap/size(deg,1) * ones(size(deg,1),1);
            obj.res0 = obj.res;
        end
        
        %% ----------------------------------------------------------------
        
        function x = get.x(obj)
            x = cat(1,obj.p(:,1));
        end
        
        function y = get.y(obj)
            y = cat(1,obj.p(:,2));
        end
        
        %% ----------------------------------------------------------------
        function l = get.l(obj)
            l = polylen(obj.p, obj.i, obj.j);
        end
        
        function lc = get.lc(obj)
            l = obj.l;
            lc = ( l + l([2:end,1]) )/2;
        end
        
        function A = get.A(obj)
            A = polyarea(obj.x, obj.y);
        end
        
        %% ----------------------------------------------------------------
        
        function dn = get.dn(obj)
            p = obj.p;                  
            i = obj.i;   % index of 3 particles on a angle
            j = obj.j;   % j = i + 1
            k = obj.k;   % k = j + 1 = i + 2
            
            xi = p(i,1); yi = p(i,2);
            xj = p(j,1); yj = p(j,2);
            xk = p(k,1); yk = p(k,2);
            
            xji = (xj + xi)/2; yji = (yj + yi)/2;
            xjk = (xj + xk)/2; yjk = (yj + yk)/2;

            dx = xji-xjk; dy = yjk-yji;
            dn = [dy, dx]./sqrt(dy.^2+dx.^2);
        end
 
        %% ----------------------------------------------------------------
        
        function force(obj, fext)
            if nargin==1; fext = 0; end
            obj.f = fext + fspring(obj.kl, obj.l0, obj.p)      + ...
                           fbending(obj.kb, obj.theta0, obj.p) + ...
                           farea(obj.ka, obj.A0, obj.p)        + ...
                           fviscous(obj.gamma, obj.v)          + ...
                           fmolecules(obj.mols, obj.dn);
        end
        
        %% ----------------------------------------------------------------
        
        function update(obj, dt, ispin)
            if nargin==2; ispin=false; end

            if ispin
                [xc0, yc0] = polyctrd(obj.p(:,1), obj.p(:,2));
            end
            
            v = obj.f/obj.mu;
            obj.p = obj.p + (obj.v + v)/2*dt;
            
            if ispin
                [xc, yc] = polyctrd(obj.p(:,1), obj.p(:,2));
                obj.p = obj.p - [xc, yc] + [xc0, yc0];
            end
            
            obj.v = v;
        end
        
        %% ----------------------------------------------------------------
       
        function cellcsk(obj)
            obj.pcsk(1:length(obj.p),:) = obj.p;  % node coordinates
            obj.pcsk = obj.Mcsk * obj.pcsk;
            [~,a1,a2,a3] = pdetrg(obj.pcsk', obj.tcsk');
            deg = acotd(-[a1', a2', a3']*2);
            obj.Acsk = mod(deg,180);
            
            if isempty(obj.Acsk0); obj.Acsk0 = obj.Acsk; end
            
            obj.dAcsk = obj.Acsk - obj.Acsk0;
            
            obj.pAcsk = obj.Mt2p * abs(obj.dAcsk(:));
        end
        
        %% ----------------------------------------------------------------
        
        function csk2mol(obj, dt)
            config('csk2mol', obj, dt);
        end
        
        %% ----------------------------------------------------------------
        function cellbdry(obj, bdry)
            if isempty(bdry); return; end
            
            dp = 1e-10;
            DistFnc = bdry.dist;
            
            % distance from p to bdrys
            d = DistFnc(obj.p);
            % Find points that are about to or have exceeded the bdrys
            I = find(d>=0);
            
            if ~isempty(I)
                % Normal and tangential direction of the local bdrys
                nx = (DistFnc(obj.p+[dp,0])-d)/dp;
                ny = (DistFnc(obj.p+[0,dp])-d)/dp;
                
                n = [nx, ny]./sqrt(nx.^2+ny.^2);
                e = [n(:,2), -n(:,1)];
                
                mu = bdry.mu(obj.p(I,:)); % Coefficient of friction
                
                % Normal force on the bdrys by point I
                N = max(sum(obj.f(I,:) .*n(I,:), 2), 0);
                
                % Tangential force to the boundary by point I
                T = sum(obj.f(I,:) .*e(I,:), 2);
                
                % friction force to point I by bdrys
                f  = min(abs(T), mu.*N).*sign(T);
                
                % Subtract tangential and normal forces from the boundary
                obj.f(I,:) = obj.f(I,:) - f.*e(I,:) - N.*n(I,:);
                obj.v(I,:) = sum(obj.v(I,:) .*e(I,:),2).*e(I,:);
                obj.p(I,:) = obj.p(I,:) - n(I,:).*d(I);
            end
        end
        
        %% ----------------------------------------------------------------
        
        function cellcell(objs)
            for i = 1:length(objs); objs(i).isadh(:) = false; end
            
            for i = 1:length(objs)
                oi = objs(i);
                for j = i+1:length(objs)
                    oj = objs(j);
                    in = inpolygon(oi.x, oi.y, oj.x, oj.y);

                    for m = find(in)'
                        r = sum((oi.p(m,:) - oj.p).^2, 2);
                        [~,n] = min(r);
                        oi.p(m,:) = (oi.p(m,:) + oj.p(n,:))/2;
                        oj.p(n,:) = oi.p(m,:);
                        oi.v(m,:) = (oi.v(m,:) + oj.v(n,:))/2;
                        oj.v(n,:) = oi.v(m,:);
                        
                        oi.isadh(m) = true; 
                        oj.isadh(n) = true;
                    end
                    
                    % adhesion
                    rx = oi.p(:,1) - oj.p(:,1)';
                    ry = oi.p(:,2) - oj.p(:,2)';
                    
                    vx = oi.v(:,1) - oj.v(:,1)';
                    vy = oi.v(:,2) - oj.v(:,2)';
                    
                    rv = rx.*vx + ry.*vy;
                    r2 = rx.*rx + ry.*ry;
                    v2 = vx.*vx + vy.*vy;
                    
                    [m,n] = find(r2<1e-3& v2<1 &rv>=0);
                    oi.p(m,:) = (oi.p(m,:) + oj.p(n,:))/2;
                    oj.p(n,:) = oi.p(m,:);
                    oi.v(m,:) = (oi.v(m,:) + oj.v(n,:))/2;
                    oj.v(n,:) = oi.v(m,:);
                end
            end
        end
        
        %% ----------------------------------------------------------------
        
        function plot(objs, cid)
        % cid = numeric (mol ID) / str (mol name) / str (field of cell)
            for j = 1:length(objs)
                obj = objs(j); 
                
                xo = obj.x; yo = obj.y;
                [xi, yi] = polyexpand(obj.x, obj.y, -0.5); % Expand a polygon
                
                x = [xo([1:end,1]); xi([1, end:-1:1])]; 
                y = [yo([1:end,1]); yi([1, end:-1:1])]; 

                xcsk0 = obj.pcsk(:,1); ycsk0 = obj.pcsk(:,2);
                ecsk = obj.ecsk;
                xcsk = [xcsk0(ecsk)'; nan(1,length(ecsk))]; 
                ycsk = [ycsk0(ecsk)'; nan(1,length(ecsk))]; 
                
                if nargin==1
                    c = zeros(size(xo));
                else
                    if isnumeric(cid)
                        c = obj.mols(cid).conc;
                    else
                        mid = find(strcmp({obj.mols.name}, cid));
                        if isempty(mid)
                            c = eval(['obj.', cid]);
                        else
                            c = obj.mols(mid).conc;
                        end
                    end
                end
                
                c = c([1:end,1,1,end:-1:1]);
                
                pAcsk = obj.pAcsk/max(obj.pAcsk)*(max(c) - min(c))+ min(c);
                ind = obj.rcsk<0.8;
                
                [xc, yc] = polyctrd(xo, yo);
                obj.xc(end+1) = xc; obj.yc(end+1) = yc;
                
                if isempty(obj.h) || ~isvalid(obj.h)
                    obj.h = patch(x, y, c, 'edgecolor','none'); hold on
                    obj.hcsk(1) = plot(xcsk(:),ycsk(:), 'color', 0.5*[1,1,1]);
                    obj.hcsk(2) = scatter(xcsk0(ind), ycsk0(ind), 5, pAcsk(ind),'filled');
                    obj.hc = plot(obj.xc, obj.yc, '-k', 'linewidth',2);
                else
                    set(obj.h, 'xdata',x, 'ydata',y, 'cdata',c)
                    set(obj.hcsk(1), 'xdata',xcsk(:), 'ydata',ycsk(:))
                    set(obj.hcsk(2), 'xdata',xcsk0(ind), 'ydata',ycsk0(ind), 'cdata', pAcsk(ind));
                    set(obj.hc,'xdata',obj.xc, 'ydata',obj.yc)
                end
                
            end
        end
        
        %% ----------------------------------------------------------------
        
    end
    
end

%% ------------------------------------------------------------------------

function [p0, l0, theta0, A0, kl, kb, ka, gamma, mu] = initcell(R, N)
% INITCELL parameters initialization of 2D cell partical model 
%

if nargin==0; N = 30; R = 3; end

beta = [0 : (2*pi/N) : (2*pi-2*pi/N)]';
p0 = R * [cos(beta), sin(beta)];

kb = 200;                     % bending constant 200
kl = 4000;                    % spring constant  2
ka = 1e6;                  % penalty coefficient of area changes

% kb = 200;                     % bending constant 200
% kl = 2000;                    % spring constant  2
% ka = 1e5;                  % penalty coefficient of area changes

gamma = 10;
mu = 200;
l0 = polylen(p0);
theta0 = polyangle(p0);
theta0(:) = 0;
A0 = polyarea(p0(:,1), p0(:,2));
end

%% ------------------------------------------------------------------------

function crossij = cross2d(ei, ej)
% CROSS2D 2D Vector Cross Product. 
%

xi = ei(:,1); yi = ei(:,2);
xj = ej(:,1); yj = ej(:,2);

crossij = xi.*yj - xj.*yi;
end

%% ------------------------------------------------------------------------

function l = polylen(p, i, j)
% POLYLEN edge lengthes of polygon
%

if nargin==1
    N = length(p);
    i = [N 1:N-1];                          % index of 2 point on a edge
    j = [1 2:N  ];                          % j = i + 1
end

l = sqrt( sum((p(j,:)- p(i,:)).^2, 2) ); 
end

%% ------------------------------------------------------------------------

function theta = polyangle(p)
% POLYANGLE angles of polygon
% http://cn.mathworks.com/matlabcentral/newsreader/view_thread/132921
%

N = length(p);

j = [1 2:N  ];                          % index of 2 point on a edge
k = [2 3:N 1];                          % k = j + 1

dx = p(k,1) - p(j,1);
dy = p(k,2) - p(j,2);

theta = atan2(dy,dx);
theta = pi - mod(pi-diff(theta([N 1:N])), 2*pi);
end

%% ------------------------------------------------------------------------

function fn = fmolecules(mols, dn)
fg = [mols.conc] * [mols.a2f]';
fn = fg.*dn;

end
%% ------------------------------------------------------------------------
function f = fspring(kl, l0, p, i, j)
% FSPRING forces on two particles connected by a spring. The two-particle 
% or spring interactions is defined as       
%
%                     V =  1/2 * kl * (l/l0 -1).
%
% The above potential energy yield the following forces
%
%           kl * (l0-l)
%     fi = ------------ * (xi - xj, yi - yj) = alpha * (xi - xj, yi - yj)
%            l0^2 * l
%
%     fj = - fi
%

if nargin==3
    N = length(p);
    i = [N 1:N-1];                          % index of 2 particles on a spring
    j = [1 2:N  ];                          % j = i + 1
end
    
f = zeros(size(p));
eji = p(j,:)- p(i, :);                  % eji = (xj-xi, yj-yi)
l = sqrt( sum((eji).^2, 2) );           % springs length

alpha = kl*(l0-l)./(l0.^2.*l);

fji = alpha .* eji;                     % fji = alpha * eij

f(i,:) = f(i,:) - fji;
f(j,:) = f(j,:) + fji;
end

%% ------------------------------------------------------------------------

function f = fbending(kb, theta0, p, i, j, k)
% FBENDING forces on triplets of particles due to angle interactions. The
% angle interaction is defined as
%
%    V = 1/2 * kb * tan(theta/2 - theta0/2)^2 = 1/2 * kb * tanhf^2
%
% The above potential energy yield the following forces
%       
%          alpha * (eji*ekj)                alpha  
%   fi =  ------------------ * eji  -   -------------- * ekj
%           |eji|^3 * |ekj|              |eji| * |ekj|
%               
%          alpha * (eji*ekj)                alpha  
%   fj =  ------------------ * ejk  -   -------------- * eij
%           |ekj|^3 * |eji|              |ekj| * |eji|
%
%   fk = -fi -fj;
%
% where alpha = 1/2 * kb * tanhalf * (1+tanhf^2) / sqrt(1-cos(theta)^2); 
% eji = -eij = (xj-xi, yj-yi);   ekj = - ejk = (xk-xj, yk-yj);
%

if nargin==3
    N = length(p);
    i = [N 1:N-1];                          % index of 3 particles on a angle
    j = [1 2:N  ];                          % j = i + 1
    k = [2 3:N 1];                          % k = j + 1 = i + 2
end

f = zeros(size(p));
ej = p(j,:)-p(i,:);                     % ej = eji = (xj-xi, yj-yi)
ek = p(k,:)-p(j,:);                     % ek = ekj = (xk-xj, yk-yj)

dotjk = sum(ej.*ek, 2);                 % ej * ek = dot(ej', ek')';
lenej = sqrt( sum(ej.^2, 2) );          % |ej| = sqrt(dot(ej',ej')');
lenek = sqrt( sum(ek.^2, 2) );          % |ek| = sqrt(dot(ek',ek')');
crojk = cross2d(ej,ek);                 % ej x ek

C0 = cos(theta0);                       
S0 = sin(theta0);
C = dotjk./(lenej.*lenek);              % cos(theta) = ej*ek /(|ej|*|ek|)
S = sign(crojk).*sqrt(1-C.^2);          % sin(theta) = sqrt(1-cos(theta))

tanhf = (S-S0) ./ (C+C0);               % tan(theta/2 - theta0/2)^2
alpha = 1/2*kb*tanhf.*(tanhf.^2 + 1)./(S+eps*(S==0));

fi =  alpha.*dotjk./(lenej.^3.*lenek) .* ej - alpha./(lenej.*lenek) .* ek;
fk = -alpha.*dotjk./(lenek.^3.*lenej) .* ek + alpha./(lenek.*lenej) .* ej;

f(i,:) = f(i,:) + fi;
f(k,:) = f(k,:) + fk;
f(j,:) = f(j,:) - fi - fk;
end

%% ------------------------------------------------------------------------

function f = farea(ka, A0, p, i, j)
% FAREA forces on two of particles due to area conservation constraint. The
% area conservation constraint is defined as
% 
%                V = 1/2 * ka * (A/A0 - 1)^2
%
% The above potential energy yield the following forces
%
%                   ka * (A - A0)
%          fi = -  --------------- * ( yj, -xj) = alpha * ( yj, -xj)
%                     2 * A0^2
%
%                   ka * (A - A0)
%          fj = -  --------------- * (-yi, -xi) = alpha * (-yi, -xi)
%                     2 * A0^2
%

if nargin==3
    N = length(p);
    i = [N 1:N-1];                          % index of 3 particles on a angle
    j = [1 2:N  ];                          % j = i + 1
end

f = zeros(size(p));

A = polyarea(p(:,1), p(:,2));           % area = sum (pi x pj)/2

alpha = -ka/2 * (A-A0)/A0^2;

fi = alpha*[ p(j,2), -p(j,1)];          % fi = alpha * ( yj, -xj)
fj = alpha*[-p(i,2),  p(i,1)];          % fj = alpha * (-yi,  xi)

f(i,:) = f(i,:) + fi;
f(j,:) = f(j,:) + fj;
end

%% ------------------------------------------------------------------------

function f = fviscous(gamma, v, i, j)
% FVISCOUS viscous force of membrane. The viscous force is defined as
%
%         fi = - gamma * (vi - vj);   fj = - fi
%

if nargin==2
    N = length(v);
    i = [N 1:N-1];                          % index of 3 particles on a angle
    j = [1 2:N  ];                          % j = i + 1
end

f = zeros(size(v));

dv = v(i,:)-v(j,:);                     % dv = vi - vj

f(i,:) = f(i,:) - gamma*dv;             % fi = - gamma * (vi - vj);
f(j,:) = f(j,:) + gamma*dv;             % fj =   gamma * (vi - vj);
end

%%-------------------------------------------------------------------------

function [x, y, xnorm, ynorm] = polyexpand(x0, y0, d)
%POLYGEXPAND Expand a polygon by a given (signed) distance

n = length(x0);
i = 1:n; j = [n, 1:n-1];

% Calculate tangent vectors
dx = x0(j) - x0(i); dy = y0(j) - y0(i);
lsq = dx.^2 + dy.^2;
dx = dx ./ max(lsq,eps);
dy = dy ./ max(lsq,eps);

Dx = zeros(n,1);    Dy = zeros(n,1);
Dx(i) = Dx(i) + dx; Dy(i) = Dy(i) + dy;
Dx(j) = Dx(j) + dx; Dy(j) = Dy(j) + dy;

% Normalize the normal
L = sqrt(Dx.^2 + Dy.^2);
xnorm = -Dy./L;   ynorm =  Dx./L;

% expand polygon (x0, y0)
x = x0 + d.*xnorm;
y = y0 + d.*ynorm;
end

%% ------------------------------------------------------------------------

function [xc, yc] = polyctrd(x, y)
N = size(x,1); i = [N 1:N-1]; j = [1 2:N  ];  % j = i + 1

xi = x(i); yi = y(i);
xj = x(j); yj = y(j);

ai = xi.*yj - xj.*yi;
area6 = 3 * sum(ai);                          % area = 1/2 * sum(ai)

xc = sum( ai.*(xi + xj) ) / area6;
yc = sum( ai.*(yi + yj) ) / area6;
end
%% ------------------------------------------------------------------------

function [p, t] = distmesh2d(fd, fh, h0, bbox, pfix)
% DISTMESH2D Mesh generation using distance functions.
% DistMesh2D is a simple surface in 2D mesh generation algorithm using 
% distance functions to define geometries.
%
%   [ P, T, STAT ] = DISTMESH2D( FD, FH, H0, BBOX, PFIX)
%
%   Input:
%
%      FD:        Distance function d(x,y,(z))
%      FH:        Scaled edge length function h(x,y,(z))
%      H0:        Initial edge length
%      BBOX:      Bounding box [xmin,ymin,(zmin); xmax,ymax,(zmax)]
%      PFIX:     Fixed node positions (N_P_FIX x 2/3)
%
%   Output:
%
%      P:         Grid vertex/node coordinates (N_P x 2/3)
%      T:         Triangle indices (N_T x 3)
%
% This function is simplified from the following code:
%   https://popersson.github.io/distmesh

% Initialization and meshing parameters.
IALG    = 2;            % Optimized algorithm selection.
IT_MIN  = 20;           % Minimum number of iterations.
IT_MINC = 50;           % Minimum number of iter. after which to call 
                        % constraint function.
                        
IT_PRT  = 25;           % Output every IT_PRT iterations.
N_RECV  = 2;            % Number of recovery iteration steps to move points
                        % outside back to boundary.
dim  = 2;

dp_tol   = -0.001*h0;   % Abs point rejection tol (dist(p)>=dp0_tol)
dtrm_tol = -0.001*h0;   % Abs dist tol for tri rejection 
                        % (t(dist(p_tcent)>=dtrm_tol) are rejected).
                        
rt_tol   =  0.3;        % Rel fraction of h0 to trigger retriangulation.
F_scale  =  1.2;        % Rel force scaling factor.
F_DCF    =  2;          % Fraction of L to L_target to allow.
dp_scale =  0.2;        % Rel fraction of computed new distance to move 
                        % points in update step.

dpc_tol = 0.001*h0;     % Abs tol for grid point movements during 
                        % convergence check.

gradeps = sqrt(eps)*h0; % Gradient computation offset.

% Initial grid point distribution, p, confined to the bounding box.
pinit{1} = bbox(1,1):h0:bbox(2,1);
pinit{2} = bbox(1,2):h0*sqrt(3)/2:bbox(2,2);
pp = cell(1,dim);
[pp{:}] = ndgrid( pinit{:} );
pp{1}(:,2:2:end) = pp{1}(:,2:2:end) + h0/2;
p = zeros(prod(size(pp{1})),dim);
for i=1:dim; p(:,i) = pp{i}(:); end

% Remove points outside the region and apply the rejection method.
p = p( call_function(fd,p)<-dp_tol, : );

if isempty(p); return; end

r0 = call_function(fh,p);   % Probability to keep point.
p  = p( rand(size(p,1),1) < min(r0)^dim./r0.^dim, : );
nfix = size(pfix,1);
if ~isempty(pfix); p = [ pfix; setdiff(p,pfix,'rows') ]; end
n_p = size(p, 1);

p0 = inf;
n_tri = 0;

while true
    % Retriangulate, if grid points have moved significantly.
    delta_p_max = max( sqrt(sum((p-p0).^2,2)) );
    if rt_tol*h0<delta_p_max
        n_tri = n_tri + 1;
        
        [p,t] = triangulate(p, fd, dtrm_tol);
        p0  = p;
        n_p = size(p,1);
        
        % Describe each edge by a unique edge_pairs of nodes.
        e = [ t(:,[1,2]); t(:,[2,3]); t(:,[3,1]) ];
        e = sort(e,2);
        e_max = max(e(:));
        if e_max*(e_max+1)<realmax
            ecomp = (e_max+1)*e(:,1) + e(:,2);
            [tmp,ind] = unique(ecomp);
            pairs = e(ind,:);
        else
            pairs = unique(e, 'rows');
        end
    end
    
    % Move mesh points based on edge lengths L and forces F.
    p1 = p(pairs(:,1),:);
    p2 = p(pairs(:,2),:);
    bars = p1 - p2;              % Bar vectors
    L = sqrt(sum(bars.^2,2));    % Bar lengths
    hbars = call_function(fh, (p1 + p2)/2);   % Rel bar mid point sizes
    % Bar target lengths.
    Ltarget = hbars*F_scale*(sum(L.^dim)/sum(hbars.^dim))^(1/dim);
    
    % Density control, remove points that are too close to each other.
    if any(Ltarget>F_DCF*L)
        p(setdiff(reshape(pairs(Ltarget>F_DCF*L,:),[],1),1:nfix),:) = [];
        n_p = size(p,1);
        p0  = inf;
        continue
    end
    
    % Compute grid point movements.
    F = max( Ltarget-L, 0 );   % Scalar bar forces.
    Fbar = F./L*ones(1,dim).*bars;
    dp = [];
    for i=1:dim
        dp = [dp, accumarray(pairs(:),[Fbar(:,i); -Fbar(:,i)],[n_p,1]) ];
    end
    
    dp(1:nfix,:) = 0;
    dp = dp_scale * dp;
    p = p + dp;
    
    % Check for convergence.
    dist = call_function( fd, p );
    delta_p_max = abs(max([sqrt(sum(dp(dist<dp_tol,:).^2,2)); -inf]));
    if delta_p_max<dpc_tol break; end
end
end

%--------------------------------------------------------------------------

function [p, t] = triangulate( p, fd, dtrm_tol )

AV_TOL = eps*1e1;   % Minimum accepted absolute area/volume.

[is_nan,tmp] = find( isnan(p) );
p(is_nan,:)  = [];

% Generate triangulation for grid points p.
t = delaunay_triangulation( p);

% Calculate simplex centers.
pc = [mean(reshape(p(t,1),size(t)),2), mean(reshape(p(t,2),size(t)),2)];

% Remove simplices with center outside region.
dist = call_function(fd,pc);
t = t(dist<dtrm_tol,:);

% Reorient simplices.
av = simpvol( p, t );
ix_flip = av<0;
t(ix_flip,[1,2]) = t(ix_flip,[2,1]);

% Remove simplices with volume < AV_TOL.
t(abs(av)<AV_TOL,:) = [];

if isempty(t); t = delaunay_triangulation(p); end
end

%--------------------------------------------------------------------------

function t = delaunay_triangulation( p)

if exist('DelaunayTri')
    t = DelaunayTri( p(:,1), p(:,2) );
else
    t = delaunay( p );
end

if isa(t,'DelaunayTri'); t = t.Triangulation; end
if isa(t,'delaunayTriangulation'); t = t.ConnectivityList; end
end

%--------------------------------------------------------------------------

function v = simpvol( p, t )
% SIMPVOL Simplex volume.
d12 = p(t(:,2),:) - p(t(:,1),:);
d13 = p(t(:,3),:) - p(t(:,1),:);
v = ( d12(:,1).*d13(:,2) - d12(:,2).*d13(:,1) )/2;
end

%--------------------------------------------------------------------------

function [ varargout ] = call_function( fun, varargin )
varargout = cell(1,max(1,nargout(fun)));
[varargout{:}] = fun( varargin{:} );
end
