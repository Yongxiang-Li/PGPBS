function [ fit ] = fit_PGPBS_Integer( X, Y, regr, corr, bs, lob, upb)
% Periodic Gaussian Process Modeling Using Bspline
    t = now;
    bsplines = bs;    Nk = [bsplines(:).p];    r = length(bsplines);
    data = struct('corr',corr, 'regr',regr, 'bs',bs, 'X',X, 'Y',Y, 'F',regr(X));
    data.beta = (data.F' * data.F) \ (data.F' * data.Y);
    
    u = linspace(0,1,5);    [u1, u2] = meshgrid(u(2:end-1));
    paras = lob(2:end) + (upb(2:end)-lob(2:end)).* [u1(:) u2(:)];
    T = unique(X(X>=lob(1) & X<=upb(1))); %T = unique([T; [(T(1:end-1)+T(2:end))/2]]);
    objk = nan(size(T,1),r); minObj = inf(r,1);
    for k = 1 : r
        data.bs = bsplines(k);
        for i = 1 : length(T)
            thetai = [repmat(T(i),size(paras,1),1) paras];
            obji = obj_func(thetai', data);
            objk(i,k) = mean(obji);
            if minObj(k)>min(obji), minObj(k)=min(obji); end
        end % figure; plot(T, -objk(:,k))
    end
    objs = mean(objk,2); % figure; plot(T, -objs)
    [~, minIndex] = min(objs);   int_period = T(minIndex);
    int_P = T;    int_likelihood = - objs;
    time = (now-t)*24*3600;
     
    [~, minIndex] = min(objs);
    AIC = 2*Nk + minObj'; % AIC = 2*Nk + objk(minIndex,:);
    [~, k] = min(AIC);   data.bs = bsplines(k);

    % exploit
    t = now;
    lob(1) = T(minIndex,1);    upb(1) = T(minIndex,1);
    theta0 = (lob+upb)/2;
    X = mod(data.X, theta0(1));    B = data.bs.func(X/theta0(1));
    gamma = ( B* B'+1e-5*var(data.Y-data.F*data.beta)*eye(data.bs.p)) ...
        \ (B*(data.Y-data.F*data.beta));
    Yhat = data.F*data.beta + B'*gamma;
    theta0(3) = std((data.Y-Yhat).^2)/std(data.Y - data.F*data.beta);
    [theta, loglik, fit, ~] = boxmin(@obj_func, theta0, lob, upb, data);
    time = time + (now-t)*24*3600;
    
    fit.theta0 = theta0;
    fit.thetahat = theta;
    fit.int_period = int_period;
    fit.int_time = time;
    fit.int_P = int_P;
    fit.int_likelihood = int_likelihood;
    fit.period = theta(1);
    fit.P = T;
    fit.likelihood = -objs;
end

function [obj, fit] = obj_func(para, data) % likelihood at interger 
    if any(para(1,:)~=para(1,1)), error('wrong input');  end
    T = para(1,1);    bs = data.bs;    
    n = size(data.Y,1);   m = bs.p;
    X = mod(data.X, T);    data.B = bs.func(X/T);    data.BtB = data.B*data.B';
    data.Yb = data.B*data.Y;    data.Fb = data.B*data.F;
    Y = data.Y - data.F*data.beta;    Yb = data.Yb - data.Fb*data.beta;
    obj = nan(size(para,2),1);   sigma2 = nan(size(para,2),1);
    for i = 1 : length(obj)
        phi = para(2,i);    delta2 = para(3,i).^2;
        r = data.corr(0, phi, bs.p, (0.5:bs.p)', 0.5);
        data.e = abs(fft(r)) + sqrt(eps*m);
        inv_R = toeplitz(sparse(delta2 * ifft(1./data.e)));
        
        varXi = (inv_R + data.BtB)/n;    L = chol(varXi, 'lower');
        YY = Y'*Y;    inv_Yb = L \ Yb;    invYY = (inv_Yb'*inv_Yb)/n;    
        sigma2(i) = (YY - invYY)/(n*delta2);
        obj(i) = m*log(n) + n*log(sigma2(i)) + (n-m)*log(delta2) + sum(log(data.e)) + 2*sum(log(diag(L)));
    end
    if nargout > 1
        gamma_E = L'\inv_Yb/n;
        gamma_V = sigma2*(eye(m) - (L' \ (L \ (data.BtB/n))))*toeplitz(r);
        fit = struct('data',data, 'period',T, 'sigma2',sigma2, ...
             'gamma_E',gamma_E, 'gamma_V',gamma_V, 'delta',para(3,:));
    end
end

function  [t, f, fit, perf] = boxmin(objfunc, t0, lo, up, data)
%BOXMIN  Minimize with positive box constraints

    % Initialize
    [t, f, fit, itdata] = start(objfunc, t0, lo, up, data);
    if  ~isinf(f)
      % Iterate
      p = length(t);
      if  p <= 2,  kmax = 2; else,  kmax = min(p,4); end
      for  k = 1 : kmax
        th = t;
        [t, f, fit, itdata] = explore(objfunc, t, f, fit, itdata, data);
        [t, f, fit, itdata] = move(objfunc, th, t, f, fit, itdata, data);
      end
    end
    perf = struct('nv',itdata.nv, 'perf',itdata.perf(:,1:itdata.nv));
end

function  [t, f, fit, itdata] = start(objfunc, t0, lo, up, data)
% Get starting point and iteration dataameters

    % Initialize
    t = t0(:);  lo = lo(:);   up = up(:);   p = length(t);
    D = 2 .^ ([1:p]'/(p+2)/2); D(1) = 1+(up(1)-lo(1))/sqrt(length(data.X))/t0(1);
    ee = find(up == lo);  % Equality constraints
    if  ~isempty(ee)
      D(ee) = ones(length(ee),1);   t(ee) = up(ee); 
    end
    ng = find(t < lo | up < t);  % Free starting values
    if  ~isempty(ng)
      t(ng) = (lo(ng) .* up(ng).^7).^(1/8);  % Starting point
    end
    ne = find(D ~= 1);

    % Check starting point and initialize performance info
    [f  fit] = objfunc(t,data);   nv = 1;
    itdata = struct('D',D, 'ne',ne, 'lo',lo, 'up',up, ...
      'perf',zeros(p+2,200*p), 'nv',1);
    itdata.perf(:,1) = [t; f; 1];
    if  isinf(f)    % Bad dataameter region
      return
    end

    if  length(ng) > 1  % Try to improve starting guess
      d0 = 16;  d1 = 2;   q = length(ng);
      th = t;   fh = f;   jdom = ng(1);  
      for  k = 1 : q
        j = ng(k);    fk = fh;  tk = th;
        DD = ones(p,1);  DD(ng) = repmat(1/d1,q,1);  DD(j) = 1/d0;
        alpha = min(log(lo(ng) ./ th(ng)) ./ log(DD(ng))) / 5;
        v = DD .^ alpha;   tk = th;
        for  rept = 1 : 4
          tt = tk .* v; 
          [ff  fitt] = objfunc(tt,data);  nv = nv+1;
          itdata.perf(:,nv) = [tt; ff; 1];
          if  ff <= fk 
            tk = tt;  fk = ff;
            if  ff <= f
              t = tt;  f = ff;  fit = fitt; jdom = j;
            end
          else
            itdata.perf(end,nv) = -1;   break
          end
        end
      end % improve

      % Update Delta  
      if  jdom > 1
        D([1 jdom]) = D([jdom 1]); 
        itdata.D = D;
      end
    end % free variables

    itdata.nv = nv;
end

function  [t, f, fit, itdata] = explore(objfunc, t, f, fit, itdata, data)
% Explore step

    nv = itdata.nv;   ne = itdata.ne;
    for  k = 1 : length(ne)
      j = ne(k);   tt = t;   DD = itdata.D(j);
      if  t(j) == itdata.up(j)
        atbd = 1;   tt(j) = t(j) / sqrt(DD);
      elseif  t(j) == itdata.lo(j)
        atbd = 1;  tt(j) = t(j) * sqrt(DD);
      else
        atbd = 0;  tt(j) = min(itdata.up(j), t(j)*DD);
      end
      [ff  fitt] = objfunc(tt,data);  nv = nv+1;
      itdata.perf(:,nv) = [tt; ff; 2];
      if  ff < f
        t = tt;  f = ff;  fit = fitt;
      else
        itdata.perf(end,nv) = -2;
        if  ~atbd  % try decrease
          tt(j) = max(itdata.lo(j), t(j)/DD);
          [ff  fitt] = objfunc(tt,data);  nv = nv+1;
          itdata.perf(:,nv) = [tt; ff; 2];
          if  ff < f
            t = tt;  f = ff;  fit = fitt;
          else
            itdata.perf(end,nv) = -2;
          end
        end
      end
    end % k

    itdata.nv = nv;
end

function  [t, f, fit, itdata] = move(objfunc, th, t, f, fit, itdata, data)
% Pattern move

    nv = itdata.nv;   ne = itdata.ne;   p = length(t);
    v = t ./ th;
    if  all(v == 1)
      itdata.D = itdata.D.^.2;
      return
    end

    % Proper move
    rept = 1;
    while  rept
      tt = min(itdata.up, max(itdata.lo, t .* v));  
      [ff  fitt] = objfunc(tt,data);  nv = nv+1;
      itdata.perf(:,nv) = [tt; ff; 3];
      if  ff < f
        t = tt;  f = ff;  fit = fitt;
        v = v .^ 2;
      else
        itdata.perf(end,nv) = -3;
        rept = 0;
      end
      if  any(tt == itdata.lo | tt == itdata.up), rept = 0; end
    end

    itdata.nv = nv;
    itdata.D = itdata.D.^.25;
end