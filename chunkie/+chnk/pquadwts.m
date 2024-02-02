function mat = pquadwts(r,d,n,d2,h,ct,bw,j,...
    rt,dt,nt,d2t,kern,opdims,t,w,opts,intp_ab,intp)
%CHNK.INTERPQUADWTS product integration for interaction of kernel on chunk 
% at targets
%
% WARNING: this routine is not designed to be user-callable and assumes 
%   a lot of precomputed values as input
%
% Syntax: mat = interpquadwts(r,d,d2,h,ct,bw,j, ...
%   rt,dt,d2t,kern,opdims,t,w,opts)
%
% Input:
%   r - chnkr nodes
%   d - chnkr derivatives at nodes
%   n - chnkr normals at nodes
%   d2 - chnkr 2nd derivatives at nodes
%   h - lengths of chunks in parameter space
%   ct - Legendre nodes at order of chunker
%   bw - barycentric interpolation weights for Legendre nodes at order of
%   chunker
%   j - chunk of interest
%   rt,dt,nt,d2t - position, derivative, normal,second derivative of select 
%               target points. if any are not used by kernel (or not well
%               defined, e.g. when not on curve), a dummy array
%               of the appropriate size should be supplied
%   kern - kernel function of form kern(srcinfo,targinfo)
%   opdims - dimensions of kernel
%   t - (Legendre) integration nodes for adaptive integration
%   w - integration nodes for adaptive integrator (t and w not necessarily 
%       same order as chunker order)
%
% Output
%   mat - integration matrix
%

if_helm2d = contains(func2str(kern),'helm2d'); % if helm2d, then kernel split quadr
if_dlp = contains(func2str(kern),'''D''') || contains(func2str(kern),'''d'''); % Hai: if dlp, is there a standard way to do strcmpi(type,'d') here?
if_slp = contains(func2str(kern),'''S''') || contains(func2str(kern),'''s'''); 

% Helsing-Ojala (interior/exterior?)
h_i = h(j);
xlohi = intp_ab*(r(1,:,j)'+1i*r(2,:,j)');         % panel end points
r_i = intp*(r(1,:,j)'+1i*r(2,:,j)');              % upsampled panel
d_i = h_i*(intp*(d(1,:,j)'+1i*d(2,:,j)'));        % r'
d2_i = h_i^2*(intp*(d2(1,:,j)'+1i*d2(2,:,j)'));   % r''
sp = abs(d_i); tang = d_i./sp;                    % speed, tangent
n_i = -1i*tang;                                   % normal
cur = -real(conj(d2_i).*n_i)./sp.^2;              % curvature
wxp_i = w.*d_i;                                     % complex speed weights (Helsing's wzp)

[mat_ho_slp, mat_ho_dlp] = SDspecialquad(struct('x',rt(1,:)' + 1i*rt(2,:)'),...
                                         struct('x',r_i,'nx',n_i,'wxp',wxp_i),xlohi(1),xlohi(2),'e');
if if_helm2d % helm2d (kernel split by J. Helsing & A. Holst: Variants of an explicit kernel-split panel-based Nystrom discretization scheme for Helmholtz boundary value problems)
  % Hai: could add helm2d slp also, but how to access naive helm2d slp & dlp matrix separately to split kernel? 
  mat_ho_slp = (-2*pi*mat_ho_slp); mat_ho_dlp = real(-2*pi*mat_ho_dlp); % Laplace slp & dlp kernel special quadr matrix 
  %mat_naive_dlp = 1i/4*(zk*abs(r_i(:).'-(rt(1,:)' + 1i*rt(2,:)'))).*besselh(1,zk*abs(r_i(:).'-(rt(1,:)' + 1i*rt(2,:)'))).*imag(wxp_i(:).'./((rt(1,:)' + 1i*rt(2,:)')-r_i(:).')); % naive helm2d dlp matrix
  targinfo = struct('r',rt); srcinfo = struct('r',[real(r_i(:))';imag(r_i(:))'],'n',[real(n_i(:))';imag(n_i(:))']);
  mat_naive = kern(srcinfo,targinfo).*abs(wxp_i)'; % naive helm2d matrix (slp or dlp or slp+dlp via kern call)
  if if_slp % slp
    mat = ( mat_naive + 2/pi*log(abs(r_i(:).'-(rt(1,:)' + 1i*rt(2,:)'))).*imag(mat_naive) ) ...
             - 2/pi*(mat_ho_slp./abs(wxp_i)').*imag(mat_naive); % Hai: kernel split formula eq(20) in J. Helsing & A. Holst
  else % dlp or slp+dlp
    mat = ( mat_naive + 2/pi*log(abs(r_i(:).'-(rt(1,:)' + 1i*rt(2,:)'))).*imag(mat_naive) + 1/(2*pi)*real(n_i(:).'./(r_i(:).'-(rt(1,:)' + 1i*rt(2,:)'))).*abs(wxp_i)' ) ...
             - 2/pi*(mat_ho_slp./abs(wxp_i)').*imag(mat_naive) - 1/(2*pi)*real(mat_ho_dlp); % Hai: kernel split formula eq(26) in J. Helsing & A. Holst, typo in Cauchy integral coefficient? should be 1/(2*pi). I think...
  end
  mat = mat*intp;
else % lap2d
  if if_dlp
    mat = real(mat_ho_dlp)*intp;  % dlp
  elseif if_slp
    mat = mat_ho_slp*intp; % slp
  else
    mat = (mat_ho_slp+real(mat_ho_dlp))*intp;  % combined layer potential
  end
end

end

function [As, Ad, A1, A2, A3, A4] = SDspecialquad(t,s,a,b,side)
% https://github.com/ahbarnett/BIE2D/blob/master/panels/LapSLP_closepanel.m
% https://github.com/ahbarnett/BIE2D/blob/master/panels/LapDLP_closepanel.m
%%% 
if nargin<5, side = 'e'; end     % interior or exterior
zsc = (b-a)/2; zmid = (b+a)/2; % rescaling factor and midpoint of src segment
y = (s.x-zmid)/zsc; x = (t.x-zmid)/zsc;  % transformed src nodes, targ pts
%figure; plot(x,'.'); hold on; plot(y,'+-'); plot([-1 1],[0 0],'ro'); % debug
N = numel(x);                            % # of targets
p = numel(s.x);                          % assume panel order is # nodes
if N*p==0
    As = 0; Ad = 0; A1=0; A2=0;
    return
end
c = (1-(-1).^(1:p))./(1:p);              % Helsing c_k, k = 1..p.
V = ones(p,p); for k=2:p, V(:,k) = V(:,k-1).*y; end  % Vandermonde mat @ nodes
P = zeros(p+1,N);      % Build P, Helsing's p_k vectorized on all targs...
d = 1.1; inr = abs(x)<=d; ifr = abs(x)>d;      % near & far treat separately
%gam = 1i;
gam = exp(1i*pi/4);  % smaller makes cut closer to panel. barnett 4/17/18
if side == 'e', gam = conj(gam); end   % note gam is a phase, rots branch cut
P(1,:) = log(gam) + log((1-x)./(gam*(-1-x)));  % init p_1 for all targs int
% upwards recurrence for near targets, faster + more acc than quadr...
% (note rotation of cut in log to -Im; so cut in x space is lower unit circle)
if N>1 || (N==1 && inr==1) % Criterion added by Bowei Wu 03/05/15 to ensure inr not empty
  for k=1:p 
    P(k+1,inr) = x(inr).'.*P(k,inr) + c(k); 
  end  % recursion for p_k
end
% fine quadr (no recurrence) for far targets (too inaccurate cf downwards)...
Nf = numel(find(ifr)); wxp = s.wxp/zsc; % rescaled complex speed weights
if Nf>0 % Criterion added by Bowei Wu 03/05/15 to ensure ifr not empty
  P(end,ifr) = sum(((wxp.*(V(:,end).*y(:)))*ones(1,Nf))./bsxfun(@minus,y,x(ifr).'));  % int y^p/(y-x)
  for ii = p:-1:2
    P( ii,ifr) = (P(ii+1,ifr)-c(ii))./x(ifr).';
  end
end
Q = zeros(p,N); % compute q_k from p_k via Helsing 2009 eqn (18)... (p even!)
% Note a rot ang appears here too...  4/17/18
%gam = exp(1i*pi/4); % 1i;  % moves a branch arc as in p_1
%if side == 'e', gam = conj(gam); end   % note gam is a phase, rots branch cut
Q(1:2:end,:) = P(2:2:end,:) - repmat(log((1-x.').*(-1-x.')),[ceil(p/2) 1]); % guessed!
% (-1)^k, k odd, note each log has branch cut in semicircle from -1 to 1
% Q(2:2:end,:) = P(3:2:end,:) - repmat(log((1-x.')./((-1-x.'))),[floor(p/2) 1]);  % same cut as for p_1
Q(2:2:end,:) = P(3:2:end,:) - repmat(log(gam) + log((1-x.')./(gam*(-1-x.'))),[floor(p/2) 1]);  % same cut as for p_1
Q = Q.*repmat(1./(1:p)',[1 N]); % k even
As = real((V.'\Q).'.*repmat((1i*s.nx)',[N 1])*zsc)/(2*pi*abs(zsc));
As = As*abs(zsc) - log(abs(zsc))/(2*pi)*repmat(abs(s.wxp)',[N 1]); % unscale, yuk
warning('off','MATLAB:nearlySingularMatrix');
% A = real((V.'\P).'*(1i/(2*pi)));         % solve for special quadr weights
Ad = ((V.'\P(1:p,:)).'*(1i/(2*pi)));         % do not take real for the eval of Stokes DLP non-laplace term. Bowei 10/19/14
%A = (P.'*inv(V))*(1i/(2*pi));   % equiv in exact arith, but not bkw stable.
if nargout>2
  R =  -(kron(ones(p,1),1./(1-x.')) + kron((-1).^(0:p-1).',1./(1+x.'))) +...
      repmat((0:p-1)',[1 N]).*[zeros(1,N); P(1:p-1,:)];  % hypersingular kernel weights of Helsing 2009 eqn (14)
  Az = (V.'\R).'*(1i/(2*pi*zsc));  % solve for targ complex-deriv mat & rescale
  A1 = Az;
  if nargout > 3
    S = -(kron(ones(p,1),1./(1-x.').^2) - kron((-1).^(0:p-1).',1./(1+x.').^2))/2 +...
        repmat((0:p-1)',[1 N]).*[zeros(1,N); R(1:p-1,:)]/2; % supersingular kernel weights
    Azz = (V.'\S).'*(1i/(2*pi*zsc^2));
    if nargout > 4
      A1 = real(Az); A2 = -imag(Az);  % note sign for y-deriv from C-deriv
      A3 = real(Azz); A4 = -imag(Azz);    
    else
      A1 = Az; A2 = Azz; 
    end
  end
end
end
