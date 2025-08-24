clear;    close all;    dbstop if error
addpath('../routines/');    rng('default')   

N = 5000;    SNR = -13;    theta = 5;
fs = 100;    T0 = 1;

[Xr ,Yr, Xg, Yg] = randPGPBS_plot(init_PBS(3,30), N, fs, T0, SNR, theta, @period_sin_gauss_cov);
Xt = Xg(1:fs*T0*2);

lob = [0.1*T0 2  0.01];
upb = [3.2*T0 10 2];

PGPBSModel = fit_PGPBS(Xr, Yr, @regpoly0, @period_sin_gauss_cov, init_BS(3,20), lob, upb);
figure(100); subplot(1,3,2);    title('(a) Regular B-spline');    hold on;
[Yt, Vt] = predict_PGPBS(PGPBSModel, Xt);
plot(Xt, Yt, '-*','Color','blue'); hold on; plot([Xt Xt], [Yt+2*sqrt(Vt) Yt-2*sqrt(Vt)], '--', 'Color','green', 'LineWidth',1)
xlabel('$t$', 'Interpreter', 'latex'); ylabel('$\hat{y}(t)$', 'Interpreter', 'latex');
legend({'Prediction Mean', '2 \times Standard Deviation'})


PGPBSModel0 = fit_PGPBS(Xr, Yr, @regpoly0, @period_sin_gauss_cov, init_PBS(3,20), lob, upb);
[Yt, Vt] = predict_PGPBS(PGPBSModel0, Xt);
subplot(1,3,3);    title('(b) Periodic B-spline');    hold on;
plot(Xt, Yt, '-*','Color','blue'); hold on; plot([Xt Xt], [Yt+2*sqrt(Vt) Yt-2*sqrt(Vt)], '--', 'Color','green', 'LineWidth',1)

xlabel('$t$', 'Interpreter', 'latex');  ylabel('$\hat{y}(t)$', 'Interpreter', 'latex');
legend({'Prediction Mean', '2 \times Standard Deviation'})
