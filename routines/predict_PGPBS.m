function [ Y, V ] = predict_PGPBS(model, X, flag)
    if nargin<3, flag = false; end
    bs = model.data.bs;
	B = bs.func(mod(X,model.period)/model.period);
    F = model.data.regr(X);   
    if flag && isfield(model, gamma_Ec)
        Y = F*model.data.beta+B'*model.gamma_Ec;
    else
        Y = F*model.data.beta+B'*model.gamma_E;
    end
    V = full(dot(B, (model.gamma_V*B)))';
end