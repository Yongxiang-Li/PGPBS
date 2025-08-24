close all;  clear;  clc;    rng('default')

addpath('..\..\routines'); 

load('outer_0.2_05_51.2k.mat')
signal = a1;
signal = signal(1:10000);
signal = signal(:);
X = (1:length(signal))'/51200;

figure; 
subplot(3,1,1); plot(X*1000, signal,'-'); set(gca,'xlim',[0 X(end)*1000]);
title('(a) Collected Vibration Signals') 
xlabel('Time (ms)') 
ylabel('Amplitude (m/s$^2$)', 'Interpreter', 'latex','FontSize', 10);

index = randperm(length(signal), length(signal)*0.9);
X = X(index);    signal = signal(index);
lob = [300/51200 1 0.1];   upb = [600/51200 40 10]; Nc = 30:5:100;

PGPBSModel = fit_PGPBS_Integer(X, signal, @regpoly1, @period_sin_gauss_cov, init_PBS(4,Nc), lob, upb);
Xt = sort(X);
[Yt, Vt] = predict_Sparse_PGPBS(PGPBSModel, Xt);

subplot(3,1,2); plot(Xt*1000, Yt,'-'); set(gca,'xlim',[0 27]);
title('(b) Reconstructed Signals of PGPBS-Int-Irr') 
xlabel('Time (ms)') 
ylabel('Amplitude (m/s$^2$)', 'Interpreter', 'latex','FontSize', 10);

PGPBSModel = fit_PGPBS(X, signal, @regpoly1, @period_sin_gauss_cov, init_PBS(4,Nc), lob, upb);
Xt = sort(X);
[Yt, Vt] = predict_Sparse_PGPBS(PGPBSModel, Xt);

subplot(3,1,3); plot(Xt*1000, Yt,'-'); set(gca,'xlim',[0 27]); hold on;
title('(c) Reconstructed Signals of PGPBS-Dec-Irr') 
xlabel('Time (ms)') 
ylabel('Amplitude (m/s$^2$)', 'Interpreter', 'latex','FontSize', 10);
