function [ Y, V ] = predict_Sparse_PGPBS(model, X)
    bs = model.data.bs;    delta2 = model.delta.^2;
    n = length(model.data.Y);    m = bs.p;
    inv_R = toeplitz(delta2 * ifft(1./model.data.e));
    rho = log(m);    v = zeros(m, 1);    u = zeros(m,1);
    varXi = inv_R + model.data.BtB + delta2*rho*eye(m);    
    L = chol(varXi, 'lower');
    threshhold = linspace(quantile(abs(model.gamma_E),0.2), quantile(abs(model.gamma_E),1),30);
    BIC = nan(size(threshhold));    df = nan(size(threshhold));    gammas = [];   
    for i = 1 : length(threshhold)
        gamma0 = zeros;    gamma = model.gamma_E;
        while max(abs(gamma-gamma0))>1e-5
            gamma0 = gamma;
            invRes = L \ (model.data.Yb-model.data.Fb*model.data.beta-delta2*(u-rho*v));
            gamma = L'\invRes;
            v = soft(gamma+u/rho, threshhold(i));
            u = u + rho*(gamma - v);
        end % figure(20); hold on; plot(gamma);
        gammas = [gammas gamma];
        df(i) = sum(v~=0); 
        BIC(i) = sum((model.data.Y-model.data.F*model.data.beta-model.data.B'*gamma).^2)...
            /model.sigma2/delta2 + log(n)*df(i);
        if df(i)==0, BIC(i) = nan; end
    end
    [~, index] = min(BIC);
    B = bs.func(mod(X, model.period) / model.period);
    F = model.data.regr(X);  
    Y = F*model.data.beta + B'*gammas(:,index);
    V = full(dot(B, (model.gamma_V*B)))';
end

function y = soft(x, T)
%   x : data (scalar or multidimensional array)
%   T : threshold (scalar or multidimensional array)
%   y : output of soft thresholding
    y = max(1 - T./abs(x), 0) .* x;
end