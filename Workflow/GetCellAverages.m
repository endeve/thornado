function [B] = GetCellAverages(Dim, nNodes, A)
% Calculates cell averages for the nodal data matrix A.
W1D   = zeros(nNodes,1);
W1    = zeros(nNodes,1);
W2    = zeros(nNodes,1);
W3    = zeros(nNodes,1);
nX    = uint64(size(A, [1, 2, 3]) / nNodes);
B     = zeros(nX);

if(nNodes == 1)
    W1D = 1.0;
elseif(nNodes == 2)
    W1D = [0.5, 0.5];
elseif(nNodes == 3)
    W1D = [5.0 / 18.0, 8.0 / 18.0, 5.0 / 18.0];
end

if(Dim == 3)
    Wq = zeros(nNodes,nNodes,nNodes);
    for i = 1:nNodes
        for j = 1:nNodes
            for k = 1:nNodes
                Wq(i,j,k) = W1D(i) * W1D(j) * W1D(k);
            end
        end
    end
    for l = 1:nX(1)
        for m = 1:nX(2)
            for n = 1:nX(3)
                B(l,m,n) = 0;
                for i = 1:nNodes
                    for j = 1:nNodes
                        for k = 1:nNodes
                            B(l,m,n) = B(l,m,n) + Wq(i,j,k) * A(i + nNodes * (l - 1), j + nNodes * (m - 1), k + nNodes * (n - 1));
                        end
                    end
                end
            end
        end
    end
elseif(Dim == 2)
    Wq = zeros(nNodes,nNodes);
    for i = 1:nNodes
        for j = 1:nNodes
            Wq(i,j) = W1D(i) * W1D(j);
        end
    end
    for l = 1:nX(1)
        for m = 1:nX(2)
            B(l,m) = 0;
            for i = 1:nNodes
                for j = 1:nNodes
                    B(l,m) = B(l,m) + Wq(i,j) * A(i + nNodes * (l - 1), j + nNodes * (m - 1));
                end
            end
        end
    end
elseif(Dim == 1)
    Wq = zeros(nNodes);
    for i = 1:nNodes
        Wq(i) = W1D(i);
    end
    for l = 1:nX(1)
        B(l) = 0;
        for i = 1:nNodes
            B(l) = B(l) + Wq(i) * A(i + nNodes * (l - 1));
        end
    end
end