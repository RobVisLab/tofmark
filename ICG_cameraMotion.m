function [K, R, C] = ICG_cameraMotion (P);
% [K, R, C] = ICG_cameraMotion (P);
% P = K * [R | -R*C];

C = ICG_normalizePoints(null(P));
[K, R] = ICG_rq(P(1:3,1:3)); 
Ksign = sign(diag(K));

sign_normalizer = diag(Ksign);
K = K*sign_normalizer;
R = sign_normalizer' * R;

K = K / K(3,3);
