function [Xr ,Yr, Xg, Yg, gamma] = randPGPBS_plot(bs, n, fs, T0, SNR, theta, corr)
    %RANDQPGP generate QPGP signal
    p = bs.p;
    r0 = corr(0, theta, p, (0.5:p)', 0.5);   
    Sr = abs(fft(r0)) + sqrt(eps*p);
    Z = randn(p, 1);
    gamma = ifft((fft(Z).*sqrt(Sr)));
    
    Xg = ((1:n)-0.5)'/fs;    Ug = bs.func(mod(Xg,T0)/T0);    Yg = Ug'*gamma;

    Xt = Xg(1:fs*T0*2); Yt = Yg(1:fs*T0*2);
    
    figure(100); subplot(1,3,1); plot(Xt, Yt);  hold on; set(gca,'ylim',[-2,5],'ytick',[-2:1:5]);
    xlabel('$t$', 'Interpreter', 'latex');  ylabel('Amplitude');  box off
    title('(a) True Signal') 
    
    noise = randn(size(Yg));
    noise = rms(Yg)/rms(noise)/10^(SNR/20) * noise;
    Yg = Yg + noise;
    Xr = ((1:n)-rand(1,n))'/fs;   Ur = bs.func(mod(Xr,T0)/T0);    Yr = Ur'*gamma;
    noise = rms(Yr)/rms(noise)/10^(SNR/20) * noise;
    Yr = Yr + noise;
end % figure; plot(Xg,Yg); hold on; plot(Xr,Yr); 

