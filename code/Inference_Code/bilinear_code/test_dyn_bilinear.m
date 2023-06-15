% Make problem size
N1 = 200;
N2 = 100;
M = 100;
K1 = 8;
K2 = 5;

% Make data
a = zeros(N1, 1);
a(randsample(N1,K1)) = 2*randn(K1, 1);
b = zeros(N2, 1);
b(randsample(N2,K2)) = 2*randn(K2, 1);


D = randn(M,N1)/sqrt(M);
F = randn(N1,N2)/sqrt(N2);
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

[a_tmp, ~, opts] = solver_L1RLS( A, y, lambda1, a_tmp, opts);

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




