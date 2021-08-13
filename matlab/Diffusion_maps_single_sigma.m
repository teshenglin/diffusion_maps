% data sampling
r = linspace(0, 1, 500)';
theta = linspace(0, 6*pi, 500)';
X = [r.*cos(theta), r.*sin(theta)];

% Visualize the original data
figure
plot(X(:,1),X(:,2),'b.','MarkerSize',12)
title('Original data')
axis equal

% Obtain the pairwise distance matrix E
% Here we use directly the square euclidean distance
E = squareform(pdist(X, 'squaredeuclidean'));
sigma = 0.1;

% Define the kernel matrix K
K = exp(-E/(sigma^2));

% Q: normalized matrix
d_sq = sqrt(sum(K,2));
Q = K./(d_sq*d_sq');

[U,S,V] = svd(Q);
V = V./(d_sq*ones(1,n));
S = diag(S);

figure
subplot(1,2,1)
plot(2:11, S(2:11), 'b.','MarkerSize',12)
title('eigenvalues 2:11')
subplot(1,2,2)
semilogy(S, 'b.')
title('eigenvalues in log')

% Project to diffusion space
Y2 = V(:,2:3).*(ones(n,1)*S(2:3)');
Y3 = V(:,2:4).*(ones(n,1)*S(2:4)');

figure
subplot(1,2,1)
plot(Y2(:,1), Y2(:,2), 'b.','MarkerSize',12)
title('DM to 2D')
axis equal
subplot(1,2,2)
plot3(Y3(:,1), Y2(:,2), Y3(:,3), 'b.','MarkerSize',12)
title('DM to 3D')
