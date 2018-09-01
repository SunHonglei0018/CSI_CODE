function n = aic_new(M, L, eigenvalues)
    aic_values = zeros(29,1);
    for ii = 1:29
        delta_n = compute_delta_n(M, ii, eigenvalues);
%         aic_values(ii) = 2*L*(M-ii)*log(delta_n)+2*ii*(2*M-ii);
%         aic_values(ii) = -(M-ii)*L*log(delta_n)+ii*(2*M-ii);
        aic_values(ii) = -(M-ii)*L*log(delta_n)+(1/2)*ii*(2*M-ii)*log(L);
    end
    [~, n] = min(aic_values); 
end

function delta_n = compute_delta_n(M, n, eigenvalues)
    sum_lam = 0;
    for ii = n+1:1:M
        sum_lam = sum_lam + eigenvalues(ii);
    end
    pro_lam = 1;
    for ii = n+1:1:M
        pro_lam = pro_lam * eigenvalues(ii);
    end
%     delta_n = ((1/(M-n))*sum_lam)/(pro_lam^(1/(M-n)));
    delta_n = (pro_lam^(1/(M-n)))/((1/(M-n))*sum_lam);
end