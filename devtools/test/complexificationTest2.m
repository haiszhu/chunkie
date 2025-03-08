%

clear
close('all')
addpath('../../')
opts = [];
opts.testfmm = false;
opts.fmmremex = false;
opts.fmmrecompile = false;
startup(opts)

zk = 2*pi;

xlims = [-10, 10];
slope = 0.5;
a_i = 3;
delta = 1;

% determine end point based on accuracy
toff = 1.5*sqrt(-log(eps))/a_i;
xc = xlims(2) + toff;
im_max = -log((eps))/real(zk)/2;
re_max = xc + max(toff, im_max/slope);

xend = re_max - xlims(2);
    
cstruct = [];
cstruct.type = 'bump';
cstruct.disc = 'adap';
maxchunklen = 2*pi/abs(zk);
cstruct.maxchunklen = maxchunklen;

[chnkrl, chnkrm, chnkrr, f] = get_boundary_curves(cstruct, xlims, ...
      delta, a_i, xend, slope);
chnkrs = [chnkrl, chnkrm, chnkrr];
chnkrtot = merge(chnkrs);

figure(1)
clf
xs = squeeze(chnkrtot.r(1,:,:)); xs = xs(:);
ys = squeeze(chnkrtot.r(2,:,:)); ys = ys(:);
plot(real(xs),ys,'k-') % what's this?
hold on;
plot(real(xs),imag(xs),'b-') % what's this?
hold on

% Solve linear system
Dk = kernel('helm', 'd', zk);

dval = -1/2;

tic, 
A = chunkermat(chnkrtot, Dk); % nystrom + correction 
A = A + eye(chnkrtot.npt)*dval;
toc


opts.corrections = true;
corr_mat = chunkermat(chnkrtot, Dk, opts); % correction


src_out = [0.1;-1]; % point source
src = [];
src.r = src_out;
Sk = kernel('helm', 's', zk); % rhs 
rhs = Sk.eval(src, chnkrtot);
maxit = 100;
tic
% tau = gmres(A, rhs, [], 1e-10, maxit);  % solve 
tau = A \ rhs;
toc

% targets & exact soln
ht = 0.1;
ytmin = min(real(chnkrtot.r(2,:)));
ytmax = ytmin + xlims(2) - xlims(1);
ytmax = ytmax/5;
ylims = [ytmin, ytmax];

xts = xlims(1):ht:xlims(2);
yts = ylims(1):ht/5:ylims(2);
[X, Y] = meshgrid(xts, yts);
[inds, ~] = inpolygon(X,Y,[real(xs);real(xs(end));real(xs(1))],...
                        [ys;[20;20]]);
targets = [X(inds).';Y(inds).'];
u_ex = Sk.eval(src, struct('r',targets));

% force smooth quadr
opts = [];
opts.forcefmm = false;
opts.forcesmooth = true;
% opts.forceadap = true;
u_tot = chunkerkerneval(chnkrtot,Dk,tau,targets,opts);

err = NaN*ones(size(X));
err(inds) = log10(abs(u_tot - u_ex));

figure(2),clf;
plot(chnkrtot, 'k.', 'LineWidth',2); hold on;
h = pcolor(X, Y, err); % shading interp
set(h, 'EdgeColor', 'none');
clim([-16,-8])
axis equal
colorbar;
set(gca,'FontSize',18);
set(gca,'FontName','CMU serif');
colormap('jet')
title('smooth quadr')

% force pquad
opts2 = [];
opts2.forcepquad = true;
opts2.side = 'i'; % left side of panel...
u_tot2 = chunkerkerneval(chnkrtot,Dk,tau,targets,opts2);

err2 = NaN*ones(size(X));
err2(inds) = log10(abs(u_tot2 - u_ex));

figure(3),clf;
plot(chnkrtot, 'k.', 'LineWidth',2); hold on;
h = pcolor(X, Y, err2); %shading interp
set(h, 'EdgeColor', 'none');
clim([-16,-8])
axis equal
colorbar;
set(gca,'FontSize',18);
set(gca,'FontName','CMU serif');
colormap('jet')
title('product quadr')

keyboard

function [chnkrl, chnkrm, chnkrr, varargout] = get_boundary_curves(cstruct, xlims, delta, ai, xend, slope)
%
%  This subroutine returns the left, right, and middle chunkers
%  for discretized curves in the complexification dirichlet paper.
%
%  Input arguments:
%    cstruct - structure containing curve info
%       cstruct.type: "bump" - wiggly gaussian bump (default)
%                     "gotham" - Fourier series that looks like gotham
%       cstruct.disc: "unif" - use chunkerfuncuni (default)
%                     "adap" - use chunkerfunc
%       cstruct.nchs: number of chunks to use on all three segments, 
%                     must be provided if using chunkerfuncuni
%                     nchs(1) will be the number of chunks on the real
%                     part, and nchs(2) will be the number of chunks
%                     on the imaginary part
%       ctstruct.maxchunklen: maximum chunk length must be provided
%                              if using chunkerfunc
%
%   xlims: limits of real part of the curve, the real part of the curve
%         is contained in (xlims(1), xlims(2))
%   delta: parameter determining support of bump function
%          bump function is supported between (xlims(1) + delta, xlims(2) - delta)
%   ai: ramp up for imaginary part which determines its centering
%   xend: additional part in complexified parts of the curve
%         the complex parts to be discretized are (xlims(2), xlims(2) +
%         xend), and (xlims(1) - xend, xlims(1))
%   slope: slope of complexified part of curve (default (1))
%
%
    type = 'bump';
    if isfield (cstruct, 'type')
        type = cstruct.type;
    end
    disc = 'unif';
    if isfield (cstruct, 'disc')
        disc = cstruct.disc;
    end

    % set up interface
    if strcmpi(type, 'bump')
        fun1 = @(t) gaus_cos(t,1,2,1,0);
    elseif strcmpi(type, 'gotham')
        as = 2*[1,-0.5,0.3];
        bs = [2*pi,pi/1.3,sqrt(2)*pi];
        cs = [0.8,1.5,3.2];

        fun1 = @(t) cos_wig(t,as,bs,cs);

    end


    % set up bump function
    a0 = 2;
    b0 = 1;
    t0 = xlims(1) + delta;
    t1 = xlims(2) - delta;
    toff = sqrt(-log(eps))/a0;
    tt0 = t0 + toff;
    tt1 = t1 - toff;
    fun2 = @(t) smooth_bump(t,a0,b0,tt0,tt1);

    % set up parameters for the imaginary part
    if nargin < 5
        slope = 1;
    end
    bi = slope/6; 
    toff = 1.5*sqrt(-log(eps))/ai;

    t0i = xlims(1) - toff;
    t1i = xlims(2) + toff;
    
    fr = @(t) prod_curve(t, fun1, fun2);

    f = @(t) prod_curve(t, fun1, fun2, ai, bi, t0i, t1i);

% Start the discretization process   
    if strcmpi(disc , 'unif')
        nchs = cstruct.nchs;
        cparams.ta = xlims(1);
        cparams.tb = xlims(2);
        chnkrm = chunkerfuncuni(fr, nchs(1), cparams);
   
        cparams.ta = xlims(1) - xend;
        cparams.tb = xlims(1);
        cparams.ifclosed = 0;
        chnkrl = chunkerfuncuni(f, nchs(2), cparams);

        cparams.ta = xlims(2);
        cparams.tb = xlims(2) + xend;
        chnkrr = chunkerfuncuni(f, nchs(2), cparams);
    elseif strcmpi(disc, 'adap')
        maxchunklen = cstruct.maxchunklen;
        cparams.ta = xlims(1);
        cparams.tb = xlims(2);
        cparams.maxchunklen = maxchunklen;
        chnkrm = chunkerfunc(fr,cparams);
   
        cparams.ta = xlims(1) - xend;
        cparams.tb = xlims(1);
        cparams.ifclosed = 0;
        chnkrl = chunkerfunc(f,cparams);

        cparams.ta = xlims(2);
        cparams.tb = xlims(2) + xend;
        chnkrr = chunkerfunc(f,cparams);
    end
    chnkrm = sort(chnkrm);
    chnkrl = sort(chnkrl);
    chnkrr = sort(chnkrr);

    varargout{1} = f;



        
end


function [f,fd,fdd] = prod_curve(t,fun1,fun2,a,b,t0,t1)

    tsz = size(t);
    t = t(:);
    f  = zeros(2,numel(t));
    fd = zeros(2,numel(t));
    fdd= zeros(2,numel(t));

    f(1,:)  = t;
    fd(1,:) = 1;
    fdd(1,:)= 0;

    [f1,f1d,f1dd] = fun1(t);
    [f2,f2d,f2dd] = fun2(t);

    f(2,:)  = f1.*f2;
    fd(2,:) = f1.*f2d + f2.*f1d;
    fdd(2,:)= f1dd.*f2 + 2*f1d.*f2d + f1.*f2dd;
 

    if (nargin < 4)
        ifcmplx = false;
    else
        ifcmplx = true;
    end
    if (ifcmplx)
    phi = @(t,u,v,z) u*(t-z).*erfc(u*(t-z))*v - exp(-u^2*(t-z).^2)/sqrt(pi)*v;
    phid= @(t,u,v,z) u*erfc(u*(t-z))*v;
    phidd=@(t,u,v,z) -u*u*exp(-u^2*(t-z).^2)*2*v/sqrt(pi);
    f(1,:) = f(1,:).' + 1i*(phi(t,a,b,t0) - phi(t,-a,b,t1)); 
    fd(1,:)= fd(1,:).' + 1i*(phid(t,a,b,t0) - phid(t,-a,b,t1));
    fdd(1,:) =fdd(1,:).' + 1i*(phidd(t,a,b,t0) - phidd(t,-a,b,t1));
    end   

    f   = reshape(f,[2,tsz]);
    fd  = reshape(fd,[2,tsz]);
    fdd = reshape(fdd,[2,tsz]);
end


function [f,fd,fdd] = gaus_cos(t,a,b,c,d)
f   = zeros(size(t));
fd  = zeros(size(t));
fdd = zeros(size(t));

f  = f  + a*exp(-t.^2/b^2).*cos(c*t+d);
fd = fd - a*exp(-t.^2/b^2).*(2*t.*cos(c*t+d)/b^2 + c*sin(c*t+d));
fdd= fdd + a*exp(-t.^2/b^2).*(4*c*t.*sin(c*t+d)/b^2 - ...
    (c^2+2/b^2-4*t.^2/b^2.*cos(c*t+d)));


end

function [f,fd,fdd] = smooth_bump(t,a,b,t0,t1)
    b = b/abs(a);
    phi=  @(t,u,v,z) u*erfc(u*(t-z))*v;
    phid= @(t,u,v,z) -u^2*exp(-u^2*(t-z).^2)*2*v/sqrt(pi);
    phidd=@(t,u,v,z)  u^4*(t-z).*exp(-u^2*(t-z).^2)*4*v/sqrt(pi);
    f  = 2-((phi(t,a,b,t0) - phi(t,-a,b,t1))); 
    fd = -((phid(t,a,b,t0) - phid(t,-a,b,t1)));
    fdd= -((phidd(t,a,b,t0) - phidd(t,-a,b,t1)));
    %f  = 0*f/2;
    %fd = 0*fd/2;
    %fdd= 0*fdd/2;
    %f = ones(size(f));
    %fd= 0*zeros(size(f));
    %fdd=0*zeros(size(f));
    
end


