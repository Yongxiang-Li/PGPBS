function v = equalspace(a, b, n)
    if n==1
        v = (a+b)/2;
    else
        v = linspace(a,b,n);
    end
end