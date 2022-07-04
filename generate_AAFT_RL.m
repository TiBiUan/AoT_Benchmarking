function Xs_final = generate_AAFT_RL (X,nsurr,PHI,GAUSS)

if (nargin<1)
    X = [];
end

[s1 s2] = size(X);

Xs_final = zeros(s1,s2,nsurr);

for j=1:nsurr
    
    % Standardise time series
    %T1 = (X-mean(X))./std(X);
    T1 = zscore(X);
    
    % Match signal distribution to Gaussian
    [Xsorted IND1] = sort(X);
    
    if ~isempty(GAUSS)
        G = GAUSS;
    else
        G = sort(randn(s1,s2));
    end
    
    T1s = zeros(size(T1));
    for i=1:s2
        T1s(IND1(:,i),i) = G(:,i);
    end
    
    % Phase randomise
    %T2 = real(ifft(abs(fft(T1s)).*exp(sqrt(-1).*angle(fft(randn(s1,s2))))));
    if ~isempty(PHI)
        T2 = CBIG_RL2017_get_PR_surrogate_RL(T1s,1,PHI);
    else
        T2 = CBIG_RL2017_get_PR_surrogate_RL(T1s,1,[]);
    end
    
    % Restore original signal distribution
    [DUM IND2] = sort(T2);
    Xs = zeros(size(X));
    for i=1:s2
        Xs(IND2(:,i),i) = Xsorted(:,i);
    end
    Xs_final(:,:,j) = Xs;
end