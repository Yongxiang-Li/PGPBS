function [ bsplines ] = init_PBS(d, Nk)
    bsplines = [];
    for i = 1 : length(Nk)
        bs.d = d;    bs.p = Nk(i);
        bs.knots = equalspace(0+0.5/bs.p,1-0.5/bs.p,bs.p)';
        bs.contrPx = equalspace(0+0.5/bs.p,1-0.5/bs.p,bs.p)';
        bs.func = @(x)PeriodicBsplineMatrix(bs.d, bs.p, bs.knots, x);
        bsplines = [bsplines; bs];
    end
end