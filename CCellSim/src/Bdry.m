classdef Bdry < matlab.mixin.Copyable
% BDRY Create Bdry (boundary) object
%
% A Bdry is a 2-D shape that can consist of one or more solid regions,
% and can have holes. Bdry objects can created from simple single geometric
% shape or combination of multiple geometric shapes
%  
% Bdry provides a function (dist) to estimate the shortest distance from 
% any point on the plane to the boundary.
%
% B = Bdry(type, mu, parms) creates a Bdry object by specifies the type, 
% friction coefficient (mu) and parameters required for the type.
% The supported types and corresponding parameters are as follows:
%
%     'line'      -  x1, x2, y1, y2       - vertex coords
%     'circle'    -  xc, yc,  r           - center and radius coords
%     'rectangle' -  x1, x2, y1, y2       - corner coords
%     'ellipse'   -  xc, yc,  a,  b, deg  - center, semi- & semi-minor 
%                                           axes and orientation
%     'polygon'   -  xv, yv               - vertex coords
%
% Methods for Boolean operations: Users can combine any of the above basic 
% shapes to form arbitrarily complex boundaries by the following operations:
%
%     -a          - Find the complement of a
%     a + b       - Find the union of a and b
%     a - b       - Find the difference of a dnd b
%     a * b       - Find the intersection of a & b, equivalent to a - (-b)
%
% Examples:
%     xv = [-35 -15 -10 10  15  35 35 15 10 -10 -15 -35]';
%     yv = [-10 -10 -2  -2 -10 -10 10 10  2   2  10  10]';
%     bdrys = Bdry('polygon'  , 0.1, xv, yv           );
%     bdrys.plot();
%
%     bdrys = Bdry('circle'   , 0.1,  0,  0, 10       ) + ...
%             Bdry('rectangle', 0.1,  0, 50, -2, 2    ) - ...
%             Bdry('circle'   , 0.1,  0,  0,  2       ) + ...
%             Bdry('circle'   , 0.1, 50,  0, 10       ) - ...
%             Bdry('rectangle', 0.0, 48, 52, -2, 2    ) + ...
%             Bdry('ellipse'  , 0.1, 25,  0,  8, 4, 0 );
%     bdrys.plot();
%
% See also POLYSHAPE, POLYBOOL.

% Zhou Lvwen: zhoulvwen@nbu.edu.cn

% October 26, 2018
    
   properties
       dist        % function (handle): return point-to-boundary distance
       mu          % function (handle): return local friction coefficient
       x           % x coordinates of boundary contour (used only for plot)
       y           % y coordinates of boundary contour (used only for plot)
       xmin        % = min(x)
       xmax        % = max(x)
       ymin        % = min(y)
       ymax        % = max(y)
   end
   
   methods
       
       function obj = Bdry(type, mu, varargin) 
       %BDRY constructor method of boundary object
           if nargin==0; return; end
           
           switch type
               case {'line', 'Line'}
                   [x1, x2, y1, y2] = deal(varargin{:});
                   dist = @(p) dLine(p, x1, x2, y1, y2);
                   xmin = min(x1,x2); xmax = max(x1,x2);
                   ymin = min(y1,y2); ymax = max(y1,y2);
                   x = [x1 x2]; y = [y1 y2];
                   
               case {'circle', 'Circle'}
                   [xc, yc, r] = deal(varargin{:});
                   dist = @(p) dCircle(p, xc, yc, r);
                   t = 0:pi/180:2*pi;
                   x = xc + r*cos(-t); y = yc + r*sin(-t);
                   xmin = xc - r; xmax = xc + r;
                   ymin = yc - r; ymax = yc + r;
                   
               case {'rectangle', 'Rectangle'}
                   [x1, x2, y1, y2] = deal(varargin{:});
                   xmin = min(x1,x2); xmax = max(x1,x2);
                   ymin = min(y1,y2); ymax = max(y1,y2);
                   dist = @(p) dRectangle(p, xmin, xmax, ymin, ymax);
                   x = [xmin, xmin, xmax, xmax, xmin];
                   y = [ymin, ymax, ymax, ymin, ymin];
                   
               case {'ellipse', 'Ellipse'}
                   [xc, yc, a, b] = deal(varargin{1:4});
                   if length(varargin)<5
                       deg = 0;
                   else
                       deg = varargin{5};
                   end
                   dist = @(p) dEllipse(p,xc,yc,a,b, deg);
                   xmin = xc-a; xmax = xc+a;
                   ymin = yc-b; ymax = yc+b;
                   S = sind(-deg); C = cosd(-deg);
                   t = 0:pi/180:2*pi;
                   x0 = a*cos(-t); y0 = b*sin(-t);
                   x = C*x0 + S*y0 + xc;
                   y = C*y0 - S*x0 + yc;
                   
               case {'polygon', 'Polygon'}
                   [x, y] = deal(varargin{:});
                   dist = @(p) dPoly(p,x,y);
                   xmin = min(x); xmax = max(x);
                   ymin = min(y); ymax = max(y);
                   
                   [x, y] = poly2cw(x, y);
                   x = x([1:end,1]); y = y([1:end,1]);
                   
               otherwise
                   error(['Only supports lines, rectangles, circles', ...
                          'ellipses, polygons. You can customize the', ...
                          'distance function (disfun) to get more', ...
                          'shape boundaries.'])
           end
           
           obj.dist = dist;
           obj.x = x;       obj.y = y;
           obj.xmin = xmin; obj.xmax = xmax;
           obj.ymin = ymin; obj.ymax = ymax;
           obj.mu = @(p) mu;
       end
       
       % ------------------------------------------------------------------
       
       function new = plus(obj, obi)
       % + PLUS union of two boundary objects
       %   C = A + B for boundary objects A and B, returns a boundary 
       %   object C that combined area of A and B.
           new = Bdry();
           new.dist = @(p) min(obj.dist(p), obi.dist(p));
           new.mu   = @(p) muplus(p, obj, obi);
           [new.x, new.y] = polybool('union',obj.x,obj.y,obi.x,obi.y);
           new.xmin = min(new.x); new.xmax = max(new.x);
           new.ymin = min(new.y); new.ymax = max(new.y);
       end
       
       % ------------------------------------------------------------------
       
       function new = minus(obj, obi)
       % - MINUS set difference of two boundary objects
       %   C = A - B for boundary objects A and B, returns a boundary 
       %   object C containing area in A but not in B.
           new = Bdry();
           new.dist = @(p) max(obj.dist(p), -obi.dist(p));
           new.mu   = @(p) muminus(p, obj, obi);
           [new.x,new.y] = polybool('subtraction',obj.x,obj.y,obi.x,obi.y);
           new.xmin = min(new.x); new.xmax = max(new.x);
           new.ymin = min(new.y); new.ymax = max(new.y);
       end
       
       % ------------------------------------------------------------------
       
       function new = uminus(obj)
       % - UMINUS complement of a boundary object
       %   B = -A for boundary objects A, returns a boundary object B
       %   containing the complement of area in B.
           new = copy(obj);
           new.dist = @(p) -obj.dist(p);
           inf = 1e10; x = [-inf -inf inf inf]; y = [-inf inf inf -inf];
           [new.x, new.y] = polybool('subtraction', x, y, obj.x, obj.y);
       end
       
       % ------------------------------------------------------------------
       
       function new = mtimes(obj, obi)
       % * MTIMES set intersection  of two boundary objects
       %   C = A * B for boundary objects A and B, returns a boundary 
       %   object C containing area common to A and B.
           new = obj - (-obi); % A * B equivalent to A - (-B)
       end
       
       % ------------------------------------------------------------------
       
       function plot(obj, method)
       %PLOT Boundary area plot.
       %   plot(B) fill the area surrounded by boudary B and draw the 
       %   boundary line
       %
       %   plot(B, method) specifies the plot method:
       %      'face'  - only fill the area surrounded by boudary B
       %      'line'  - only draw the boundary line

           if isempty(obj); return; end
           
           if nargin==1; method = 'otherwise'; end
           switch method
               case 'face'
                   [f, v] = poly2fv(obj.x, obj.y);
                    patch('Faces', f, 'Vertices', v, 'FaceColor', ...
                           [1,0.8,0.8], 'EdgeColor', 'none')
               case 'line'
                   plot(obj.x, obj.y, '-k')
               otherwise
                   [f, v] = poly2fv(obj.x, obj.y); hold on
                   patch('Faces', f, 'Vertices', v, 'FaceColor', ...
                           [1,0.8,0.8], 'EdgeColor', 'none')
                   plot(obj.x, obj.y, '-k')
           end
           axis image;
           
           xlim = get(gca, 'XLim'); ylim = get(gca, 'YLim'); 
           inf = 1e10;
           if xlim(1)<-inf/2; xlim(1) =  inf; end
           if xlim(2)> inf/2; xlim(2) = -inf; end
           if ylim(1)<-inf/2; ylim(1) =  inf; end
           if ylim(2)> inf/2; ylim(2) = -inf; end
           xmin = min(xlim(1), obj.xmin); xmax = max(xlim(2), obj.xmax);
           ymin = min(ylim(1), obj.ymin); ymax = max(ylim(2), obj.ymax);
           axis(1.25*[xmin, xmax, ymin, ymax]); 
       end
   end
    
end

% -------------------------------------------------------------------------

function mu = muminus(p, obj, obi)
% MUMINUS return local friction coefficient for set difference of two Bdrys
flag = obj.dist(p)>-obi.dist(p);
mu = obj.mu(p).*flag + obi.mu(p).*(~flag);  
end

function mu = muplus(p, obj, obi)
% MUPLUS return local friction coefficient for union of two Bdrys
flag = obj.dist(p)< obi.dist(p);
mu = obj.mu(p).*flag + obi.mu(p).*(~flag);  
end

% -------------------------------------------------------------------------

function d = dLine(p,x1,x2,y1,y2)
% DLINE Find minimum distances from points to a line. 
% By convention, a point located at the left hand side of the line
% is inside the region and it is assigned a negative distance value.
t = [y2-y1; x1-x2]; t = t./sqrt(sum(t.^2));
d = p - [x1 y1];
d = d * t;
end

% -------------------------------------------------------------------------

function d = dCircle(p, xc,yc, r)
% DCIRCLE Find minimum distances from points to a circle
d = sqrt( sum((p - [xc, yc]).^2,2) ) - r;
end

% -------------------------------------------------------------------------

function d = dRectangle(p, x1, x2, y1, y2)
% DRECTANGLE Find minimum distances from points to a rectangle
d = [ [x1,y1]-p, p-[x2,y2] ];
d = max(d,[],2);
end

% -------------------------------------------------------------------------

function d = dPoly(p, xv, yv)
% DPOLY Find minimum distances from points to a polygon.
% This function is simplified and optimized from the following code:
%   https://www.mathworks.com/matlabcentral/fileexchange/12744-distance-
%   from-points-to-polyline-or-polygon

xp = p(:,1); yp = p(:,2);

% number of points and number of vertices in polyline
nv = length(xv);  np = length(xp);

% matrix of distances between all points and all vertices
% dpv(j,k) - distance from point j to vertex k  
dpv = sqrt((xv'-xp).^2 + (yv'-yp).^2);

% Find the vector of minimum distances to vertices. 
[dpv_min, I_dpv_min] = min(abs(dpv),[],2);

% vector of distances between each pair of consecutive vertices
dx = xv([2:end,1]) - xv;
dy = yv([2:end,1]) - yv;
vds = sqrt(dx.^2 + dy.^2);

% Build the rotation matrix from original to rotated system
ctheta = dx./vds; stheta = dy./vds;
C1 = [ ctheta,  stheta]';
C2 = [-stheta,  ctheta]';

% translation
r  = p*C1 - sum([xv, yv]'.*C1);
cr = p*C2 - sum([xv, yv]'.*C2);

% Find the projections that fall inside the segments
is_in_seg = (r>0) & (r<vds');

% find the minimum distances from points to their projections that fall
% inside the segments (note, that for some points these might not exist,
% which means that closest points are vertices)
B = NaN(np,nv);
B(is_in_seg) = cr(is_in_seg);
[cr_min, I_cr_min] = min(abs(B),[],2);

% Build the logical index which is true for members that are vertices,
% false otherwise. Case of NaN in cr_min means that projected point falls
% outside of a segment, so in such case the closest point is vertex.

% point's projection falls outside of ALL segments
cond1 = isnan(cr_min); 

% point's projection falls inside some segments, find segments for which
% the minimum distance to vertex is smaller than to projection
cond2 = (I_cr_min ~= I_dpv_min) & (cr_min > dpv_min);

is_vertex = (cond1 | cond2);

% build the minimum distances vector
d = cr_min;
d(is_vertex) = dpv_min(is_vertex);

% mimic the functionality of ver. 1.0 - make all distances negative for
% points INSIDE the polygon
in = inpolygon(xp, yp, xv, yv);
d(in) = -d(in);
end

% -------------------------------------------------------------------------

function d = dEllipse(p, xc, yc, a, b, deg)
% DELLIPSE Find minimum distances from points to a ellipse.
% This function is simplified and optimized from the following code:
%   https://www.mathworks.com/matlabcentral/fileexchange/27708-distance-
%   from-points-to-an-ellipse

c = [xc, yc]; r = [ a, b];
tol = 1e-9; % tolerance

%  First handling the circle case
if abs((a-b)/a)<tol
    d = sqrt(sum((p - c).^2,2)) - mean(r); return;
end

n = size(p,1);
proj = zeros(n,2);

% Now dealing with proper ellipses
aa = a^2;  bb = b^2;
tol_a  = tol*a;
tol_b  = tol*b;
tol_aa = tol*aa;

% Matrix Q for rotating the points and the ellipse to the canonical system
S = sind(deg); C = cosd(deg);
Q = [C -S; S C];

% data points in canonical coordinates
p0  = (p - c)*Q;
p0abs = abs(p0);
Tini = max(r.* (p0abs - r), [], 2);

% main loop over the data points
for i = 1:n
    pj = p0abs(i,:);
    z = sign(pj) + (pj==0);
    u = p0abs(i,1);  v = p0abs(i,2);
    ua = u*a;       vb = v*b;

    %  does the point lie on the minor axis?
    if u<tol_a
        if p0(i,2)<0; proj(i,:) = [0, -b]; else; proj(i,:) = [0, b]; end
        continue;
    end

    % does the point lie on the major axis?
    if v<tol_b
        if u<a-bb/a
            xproj = aa*u/(aa-bb);
            proj(i,:) = z.*[xproj, b*sqrt(1-min((xproj/a)^2, 1))];
        else
            proj(i,:) = z.*[a 0];
        end
        continue;
    end

    % generic case: start the iterative procedure
    T = Tini(i);
    for iter = 1:100
        Taa = T + aa;     Tbb = T + bb;
        PP1 = (ua/Taa)^2; PP2 = (vb/Tbb)^2;
        F  = PP1 + PP2 - 1;
        if F<0; break; end
        Fder = 2*(PP1/Taa + PP2/Tbb);
        Ratio = F/Fder;
        if (Ratio<tol_aa); break; end
        T = T + Ratio;
    end
    % compute the projection of the point onto the ellipse
    xproj = p0(i,1)*aa/Taa;
    proj(i,:) = [xproj sign(p0(i,2))*b*sqrt(1-min((xproj/a)^2,1))];
end

% rotate back to the original system
proj = proj*Q' + c;
in = sign(sum((p0./r).^2,2)-1);
d = sqrt(sum((p - proj).^2,2)) .* in;
end