using LinearAlgebra
using MAT
using Plots
using Statistics

file = matread("Data2.mat");
X = file["data"];

# Visualize the first two dimension of the data
scatter(X[:,1], X[:,2], aspect_ratio=:equal, leg=false)

# input: X - data points
# output: E - distance matrix
#
function get_E(X)
    n = size(X,1); # number of points
    E = zeros(n,n); # the diagonal entries are 0

    for ii=1:n
        for jj=ii+1:n
            E[ii,jj] = norm(X[ii,:] - X[jj,:]);
            E[jj,ii] = E[ii,jj];
        end
    end

    return E;
end

# get the distance matrix D
E = get_E(X);

# evaluate sigma
E_sort = sort(E, dims=2);

k = 7;
sigma_loc = E_sort[:, k+1];

# input 1: E - distance matrix
# input 2: sigma - constant
#
# output: K - kernal matrix
#
function get_K(E, sigma_loc)
    n = size(E,1);
    K = ones(n,n);

    for ii = 1:n
        for jj = ii+1:n
            K[ii,jj] = exp(-E[ii,jj]^2/(sigma_loc[ii]*sigma_loc[jj]));
            K[jj,ii] = K[ii,jj]
        end
    end

    return K;
end

# get the kernal matrix K
K = get_K(E, sigma_loc);

# input: K - kernal matrix
#
# output 1: Q
# output 2: d_sq - sqrt{D}
#
function get_Q(K)
    n = size(K,1);
    Q = zeros(n,n);
    d_sq = zeros(n);

    # d_sq = sqrt{D}
    for ii = 1:n
        d_sq[ii] = sqrt(sum(K[ii,:]));
    end

    # get components of Q
    for ii = 1:n
        for jj = 1:n
            Q[ii,jj] = K[ii,jj]/(d_sq[ii]*d_sq[jj]);
        end
    end

    return Q, d_sq;
end

# get Q and d_sq
Q, d_sq = get_Q(K);

# input 1: Q
# input 2: d_sq - sqrt{D}
#
# output 1: v - eigenvectors
# output 2: s - eigenvalues
#
function get_eig(Q, d_sq)

    n = size(Q, 1);

    U,S,V = svd(Q); # U and S contains eigenvectors and eigenvalues repectively
                    # which is arranged in descending power.

    v = zeros(n,n);

     for ii = 1 : n
        V[ii,:] = V[ii,:]/d_sq[ii];
    end


    for ii = 1 : n
        v[:,ii] = V[:,ii]/norm(V[:,ii]);
    end

    return v, S;
end

c = 3 ; # the desired reduced dimension

v, s = get_eig(Q , d_sq);
p1 = scatter(s[1:10], label="eigenvalues 1:10");
p2 = plot(log.(s), label="eigenvalues in log");
plot(p1, p2, layout=2)

function get_Y(v, S, c)

    n = size(v,1);
    Y = zeros(n,c);

    # get components of diffusion map Y
    for ii = 1:c
        Y[:,ii] = v[:,ii+1].*S[ii+1];
    end

    return Y ;
end

Y = get_Y(v, s, c);

# print diffution map
p1 = scatter(Y[:,1], Y[:,2], label="2D", aspect_ratio=:equal)
p2 = scatter(Y[:,1], Y[:,2], Y[:,3], label="3D", aspect_ratio=:equal)
plot(p1, p2, layout=2)
