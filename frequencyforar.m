function Fstats= frequencyforar(ARmodel, kx, ky, w, cov_xy, fs)
% get the coefs for AR
coefs = ARmodel.A; % the first coef is always eye matrix; the secod for the first lag, the third for the second lag, and so on
% form the coefs matrix A
% get the coefs of A, B, C, D
% orders = [p * ones(kx,kx), p1 * ones(kx,ky); q1 * ones(ky,kx), q * ones(ky,ky)];
AR = max(max(ARmodel.na));
A = zeros(kx,kx,AR); B = zeros(kx,ky,AR);
C = zeros(ky,kx,AR); D = zeros(ky,ky,AR);
dx = kx; dy = ky;
for i = 1 : AR
    A(:,:,i) = coefs(1:kx,1:kx,i+1);
    B(:,:,i) = coefs(1:kx, kx+1:end, i+1);
    C(:,:,i) = coefs(kx+1:end, 1:kx, i+1);
    D(:,:,i) = coefs(kx+1:end, kx+1:end, i+1);
end
% fourier transform
a1x = eye(kx); a2x = zeros(kx,ky);
a1y = zeros(ky,kx); a2y = eye(ky);
for i = 1 : AR
    a1x = a1x + A(:,:,i) * exp(-1i*2*pi*w/fs*i);
    a2x = a2x + B(:,:,i) * exp(-1i*2*pi*w/fs*i);
    a1y = a1y + C(:,:,i) * exp(-1i*2*pi*w/fs*i);
    a2y = a2y + D(:,:,i) * exp(-1i*2*pi*w/fs*i);
end
inva2y = eye(ky) / a2y;
% inverse of the coefs matrix A
Hxx = eye(kx)  / ( a1x - a2x * inva2y * a1y); 
Hxy = - Hxx *  a2x * inva2y ;
Hyx = - inva2y * a1y * Hxx;
Hyy = ( inva2y * a1y * Hxx * a2x + eye(dy)) * inva2y;

H = [Hxx, Hxy; Hyx, Hyy];



% covariance matrix
sigma33 = cov_xy(1:dx,1:dx);
sigma34 = cov_xy(1:dx, dx+1:end);
sigma44 = cov_xy(dx+1:end, dx+1:end);

% %% for x -?-> y
% Px2y = [ eye(dx), -sigma34 / sigma44; zeros(dy,dx), eye(dy)];
% Hx2y = H / Px2y;
% Syy1 = Hx2y(dx+1:end, dx+1:end) * sigma44 * Hx2y(dx+1:end, dx+1:end)' /fs;
% Syy = Syy1 + Hx2y(dx+1:end, 1:dx) *....
%     (sigma33 - sigma34 / sigma44 * sigma34') *...
%     Hx2y(dx+1:end, 1:dx)' /fs;
% % Fstats(1) = log( norm(Syy, 2) / norm(Syy1, 2) );
% Fstats(1) = log( abs(trace(Syy) / trace(Syy1)) );
% % normsyy = norm(Syy,2);
% % normsyy = det(Syy);

%% for y -?-> x
Py2x = [eye(dx), zeros(dx,dy); -sigma34' / sigma33, eye(dy)];
Hy2x = H / Py2x;
Sxx1 = Hy2x(1:dx, 1:dx) * sigma33 * Hy2x(1:dx, 1:dx)' /fs;
Sxx = Sxx1 + Hy2x(1:dx, dx+1:end) * ...
    (sigma44 - sigma34' / sigma33 * sigma34 ) * ...
    Hy2x(1:dx, dx+1:end)' /fs;
% Fstats(2) = log( norm(Sxx, 2) / norm(Sxx1, 2) );
Fstats = log( abs(trace(Sxx) / trace(Sxx1)) );
% normsxx = det(Sxx);
% H1 = abs(trace(Sxx));
% H2 = abs(trace(Sxx1));






