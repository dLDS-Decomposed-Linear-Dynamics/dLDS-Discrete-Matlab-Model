% Make problem size
N1 = 100;
N2 = 200;
M = 300;
K1 = 5;
K2 = 8;

% Make data
a = zeros(N1, 1);
a(randsample(N1,K1)) = 2*randn(K1, 1);
b = zeros(N2, 1);
b(randsample(N2,K2)) = 2*randn(K2, 1);
W0 = a*b';

% Made 'dictionary' and 'measurements'
D = randn(M,N1*N2);
D = D*diag(1./sqrt(sum(D.^2,1)));
x_clean = D*W0(:);
x = x_clean + 0.0001*randn(M,1);

D2 = randn(M,N1,N2)/sqrt(M);
y_clean = sum(D2.*repmat(reshape(b,[1,1,N2]),[M,N1,1]),3)*a;
y = y_clean + 0.001*randn(M,1);

% Set up TFOCS
opts.tol = 1e-4;
opts.printEvery = 0;
lambda1 = 0.2;
lambda2 = 0.2;

% Lifting method
Af = @(z) D*vec(z);
Ab = @(z) reshape((D')*z, N1, N2);
% A = linop_handles({[M,1], [N1,N2]}, Af, Ab, 'R2R');
A = linop_handles({[N1,N2], [M,1]}, Af, Ab, 'R2R');

Wstart = zeros(N1, N2);
[W_finNS, out, opts ] = tfocs( smooth_quad, { A, -x }, prox_nuclear(1), Wstart, opts );
[W_fin, out, opts ] = tfocs( smooth_quad, { A, -x }, prox_spnuclear(0.01*[1,0.02,0.02],0), Wstart, opts );
[U_NS,S_NS,V_NS] = svds(W_finNS,1);
[U,S,V] = svds(W_fin,1);
% Plotting
% subplot(2,1,1), stem([a,U(:)*(a'*U(:))])
% subplot(2,1,2), stem([b,V(:)*(b'*V(:))])
% 
% subplot(2,1,1), stem([a/norm(a),U(:)])
% subplot(2,1,2), stem([b/norm(b),V(:)])

a_x = 1:numel(a);
b_x = 1:numel(b);

subplot(2,2,1), hold off
stem(a_x(a~=0),a(a~=0),'ob','LineWidth',3)
subplot(2,2,1), hold on
stem(a_x(U~=0),U(U~=0)*(a'*U(:)),'sr','MarkerSize',5)
title('Estimate of a', 'FontSize', 20)
set(gca,'FontSize',18)
subplot(2,2,1), hold off
subplot(2,2,3), hold off
stem(b_x(b~=0),b(b~=0),'ob','LineWidth',3)
subplot(2,2,3), hold on
stem(b_x(V~=0),V(V~=0)*(b'*V(:)),'sr','MarkerSize',5)
title('Estimate of b', 'FontSize', 20)
set(gca,'FontSize',18)
subplot(2,2,3), hold off


subplot(2,2,2), hold off
stem(a_x(a~=0),a(a~=0),'ob','LineWidth',3)
subplot(2,2,2), hold on
stem(a_x(U_NS~=0),U_NS(U_NS~=0)*(a'*U_NS(:)),'sr','MarkerSize',5)
title('Estimate of a', 'FontSize', 20)
set(gca,'FontSize',18)
subplot(2,2,2), hold off
subplot(2,2,4), hold off
stem(b_x(b~=0),b(b~=0),'ob','LineWidth',3)
subplot(2,2,4), hold on
stem(b_x(V_NS~=0),V_NS(V_NS~=0)*(b'*V_NS(:)),'sr','MarkerSize',5)
title('Estimate of b', 'FontSize', 20)
set(gca,'FontSize',18)
subplot(2,2,4), hold off


% subplot(2,2,2), stem([a,U_NS(:)*(a'*U_NS(:))])
% title('Estimate of a', 'FontSize', 20)
% subplot(2,2,4), stem([b,V_NS(:)*(b'*V_NS(:))])
% title('Estimate of a', 'FontSize', 20)


%%
% Alternating directions method
a_tmp = zeros(N1,1);
b_tmp = rand(N2,1);%/sqrt(N2);
for kk = 1:300
    a_tmp_old = a_tmp;
    b_tmp_old = b_tmp;
    
    D_tmp = sum(D2.*repmat(reshape(b_tmp,[1,1,N2]),[M,N1,1]),3);
    Af = @(z) D_tmp*z;
    Ab = @(z) (D_tmp')*z;
    A = linop_handles([M,N1], Af, Ab, 'R2R');
    
    [a_tmp, ~, opts] = solver_L1RLS( A, y, lambda1, a_tmp, opts);
    
    D_tmp = reshape(sum(D2.*repmat(reshape(a_tmp,[1,N1,1]),[M,1,N2]),2),M,N2);
    Af = @(z) D_tmp*z;
    Ab = @(z) (D_tmp')*z;
    A = linop_handles([M,N2], Af, Ab, 'R2R');
    
    [b_tmp, ~, opts] = solver_L1RLS( A, y, lambda2, b_tmp, opts);
    fprintf('Iteration = %3.0d:diff a: %d, diff b: %d, MSE: %d, Energy: %d\n', ...
        kk, norm(a_tmp-a_tmp_old)/norm(a_tmp_old), ...
        norm(b_tmp-b_tmp_old)/norm(b_tmp_old), ...
        norm(y_clean - Af(b_tmp))^2/norm(y_clean)^2, ...
        norm(y - Af(b_tmp))^2 + lambda1*sum(abs(a_tmp)) + lambda2*sum(abs(b_tmp)))
end

% Plotting

subplot(2,1,1), stem([a,a_tmp])
subplot(2,1,2), stem([b,b_tmp])



