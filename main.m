clc
clear all

fC = 2e9;
v = 30/3.6;
Ns = 10000;
Ts = 0.1e-3;
fs = 1/Ts;
fD = v/(3e8)*fC;
f = (-60:fs/Ns:60);

% Rayleigh PSD(theoretical)
S = 1./(pi*fD*(sqrt(1-(f/fD).^2)));
S(abs(f)>=fD) = 0;

% Rayleigh PSD(spectral method)
c_s = spectrummethod(fD,Ts,Ns,0);
[S_rayleigh_psd,f1] = pwelch(c_s,[],[],[],fs,'centered');
S_rayleigh_psd = S_rayleigh_psd /(sum(S_rayleigh_psd)*(fs/Ns));

% Rayleigh PSD(filter method)
c_f = filtermethod(fD,Ts,Ns,0);
[F_rayleigh_psd,f2] = pwelch(c_f,[],[],[],fs,'centered');
F_rayleigh_psd = F_rayleigh_psd /(sum(F_rayleigh_psd)*(fs/Ns));

figure(1)
plot(f,pow2db(S)),xlim([-100,100]);
hold on;
plot(f1,pow2db(S_rayleigh_psd));
hold on;
plot(f2,pow2db(F_rayleigh_psd));
title('PSD of Rayleigh theoretical, spectrum & filter method');
xlabel('Frequency (Hz)');
ylabel('PSD (dB/Hz)');
legend('theoretical','spectrum method','filter method');

% Rician PSD(theoretical,filter & spectrum method)
Kc1 = 0;
T_rician_psd1 = makedist('Rician','s',Kc1,'sigma',1);
F_rician_psd1 = fitdist(abs(filtermethod(fD,Ts,Ns,Kc1)),'Rician');
S_rician_psd1 = fitdist(abs(spectrummethod(fD,Ts,Ns,Kc1)'),'Rician');

Kc2 = 1;
T_rician_psd2 = makedist('Rician','s',Kc2,'sigma',1);
F_rician_psd2 = fitdist(abs(filtermethod(fD,Ts,Ns,Kc2)),'Rician');
S_rician_psd2 = fitdist(abs(spectrummethod(fD,Ts,Ns,Kc2)'),'Rician');

Kc3 = 10;
T_rician_psd3 = makedist('Rician','s',Kc3,'sigma',1);
F_rician_psd3 = fitdist(abs(filtermethod(fD,Ts,Ns,Kc3)),'Rician');
S_rician_psd3 = fitdist(abs(spectrummethod(fD,Ts,Ns,Kc3)'),'Rician');

figure(2)
subplot(3,1,1)
plot(T_rician_psd1);
hold on;
plot(F_rician_psd1);
hold on;
plot(S_rician_psd1);
subtitle('Kc=0')
ylabel('PSD(dB/Hz)');
title('PSD of Racian theoretical, filter & spectrum method');
subplot(3,1,2)
plot(T_rician_psd2);
hold on;
plot(F_rician_psd2);
hold on;
plot(S_rician_psd2);
subtitle('Kc=1')
ylabel('PSD(dB/Hz)');
subplot(3,1,3)
plot(T_rician_psd3);
hold on;
plot(F_rician_psd3);
hold on;
plot(S_rician_psd3);
subtitle('Kc=10')
xlabel('Frequency(Hz)');
ylabel('PSD(dB/Hz)');
legend('theoretical','','filter method','','spectrum method');

% pdf(thoretical)
pdf_t_kc1 = pdf(T_rician_psd1,f);
pdf_t_kc2 = pdf(T_rician_psd2,f);
pdf_t_kc3 = pdf(T_rician_psd3,f);

% pdf(filter method)
pdf_f_kc1 = pdf(F_rician_psd1,f);
pdf_f_kc2 = pdf(F_rician_psd2,f);
pdf_f_kc3 = pdf(F_rician_psd3,f);

% pdf(spectrum method)
pdf_s_kc1 = pdf(S_rician_psd1,f);
pdf_s_kc2 = pdf(S_rician_psd2,f);
pdf_s_kc3 = pdf(S_rician_psd3,f);

figure(3)
subplot(3,1,1)
plot(f,pdf_t_kc1),xlim([-5,5]);
hold on;
plot(f,pdf_f_kc1),xlim([-5,5]);
hold on;
plot(f,pdf_s_kc1),xlim([-5,5]);
subtitle('Kc = 0');
title('PDF of theoretical, filter & spectrum method');
subplot(3,1,2)
plot(f,pdf_t_kc2),xlim([-5,5]);
hold on;
plot(f,pdf_f_kc2),xlim([-5,5]);
hold on;
plot(f,pdf_s_kc2);
subtitle('Kc = 1'),xlim([-5,5]);
ylabel('PDF');
subplot(3,1,3)
plot(f,pdf_t_kc3),xlim([5,15]);
hold on;
plot(f,pdf_f_kc3),xlim([5,15]);
hold on;
plot(f,pdf_s_kc3),xlim([5,15]);
subtitle('Kc = 10');
xlabel('Frequency(Hz)');
legend('theoretical','filter method','spectrum method');

% cdf(thoretical)
cdf_t_kc1 = cdf(T_rician_psd1,f);
cdf_t_kc2 = cdf(T_rician_psd2,f);
cdf_t_kc3 = cdf(T_rician_psd3,f);

% cdf(filter method)
cdf_f_kc1 = cdf(F_rician_psd1,f);
cdf_f_kc2 = cdf(F_rician_psd2,f);
cdf_f_kc3 = cdf(F_rician_psd3,f);

% cdf(spectrum method)
cdf_s_kc1 = cdf(S_rician_psd1,f);
cdf_s_kc2 = cdf(S_rician_psd2,f);
cdf_s_kc3 = cdf(S_rician_psd3,f);

figure(4)
subplot(3,1,1)
plot(f,cdf_t_kc1),xlim([-10,10]);
hold on;
plot(f,cdf_f_kc1),xlim([-10,10]);
hold on;
plot(f,cdf_s_kc1),xlim([-10,10]);
subtitle('Kc = 0');
title('CDF of theoretical, filter & spectrum method');
subplot(3,1,2)
plot(f,cdf_t_kc2),xlim([-10,10]);
hold on;
plot(f,cdf_f_kc2),xlim([-10,10]);
hold on;
plot(f,cdf_s_kc2);
subtitle('Kc = 1'),xlim([-10,10]);
ylabel('CDF');
subplot(3,1,3)
plot(f,cdf_t_kc3),xlim([-0,20]);
hold on;
plot(f,cdf_f_kc3),xlim([-0,20]);
hold on;
plot(f,cdf_s_kc3),xlim([-0,20]);
subtitle('Kc = 10');
xlabel('Frequency(Hz)');
legend('theoretical','filter method','spectrum method');

% auto-correlation
%for Kc = 0
Kc_1 = 0;
t = (-Ns*Ts:Ts:Ns*Ts);
ac_t1 = (Kc_1^2) + besselj(0,2*pi*fD*t);
[ac_f1,lags_f1] = xcorr(filtermethod(fD,Ts,Ns,Kc_1),'unbiased');
[ac_s1,lags_s1] = xcorr(spectrummethod(fD,Ts,Ns,Kc_1),'unbiased');

%for Kc = 1
Kc_2 = 1;
t = (-Ns*Ts:Ts:Ns*Ts);
ac_t2 = (Kc_2^2) + besselj(0,2*pi*fD*t);
[ac_f2,lags_f2] = xcorr(filtermethod(fD,Ts,Ns,Kc_2),'unbiased');
[ac_s2,lags_s2] = xcorr(spectrummethod(fD,Ts,Ns,Kc_2),'unbiased');

%for Kc = 10
Kc_3 = 10;
t = (-Ns*Ts:Ts:Ns*Ts);
ac_t3 = (Kc_3^2) + besselj(0,2*pi*fD*t);
[ac_f3,lags_f3] = xcorr(filtermethod(fD,Ts,Ns,Kc_3),'unbiased');
[ac_s3,lags_s3] = xcorr(spectrummethod(fD,Ts,Ns,Kc_3),'unbiased');

figure(5)
subplot(3,1,1)
plot(t,ac_t1),xlim([-1,1]);
hold on;
plot(lags_s1.*Ts,ac_s1),xlim([-1,1]);
hold on;
plot(lags_f1.*Ts,ac_f1),xlim([-1,1]);
subtitle('Kc = 0');
title('ACF of theoretical, filter & spectrum method');
subplot(3,1,2)
plot(t,ac_t2),xlim([-1,1]);
hold on;
plot(lags_s2.*Ts,ac_s2),xlim([-1,1]);
hold on;
plot(lags_f2.*Ts,ac_f2),xlim([-1,1]);
subtitle('Kc = 1');
ylabel('Autocorrelation');
subplot(3,1,3)
plot(t,ac_t3),xlim([-1,1]);
hold on;
plot(lags_s3.*Ts,ac_s3),xlim([-1,1]);
hold on;
plot(lags_f3.*Ts,ac_f3),xlim([-1,1]);
subtitle('Kc = 10');
xlabel('Lag');
legend('theoretical','spectrum method','filter method');

% Time and Frequency-Varying Rician Fading Channels
N = 300;
M = 64;
fdTs = [0.1 0.005];
Kc_value = [0 1 10];

for L = 1:3
    for i = 1:3
        for j = 1:2
            for k = 1: L
                c_spec(k,:) = spectrummethod(fdTs(j),Ts,N,Kc_value(i));
            end    
            [l_spec, ~] = size(c_spec);
            c_padded_spec = [c_spec; zeros(M-l_spec,N)];
            C_S = fft(c_padded_spec,M);
            figure()
            subplot(1,2,1)
            mesh(0:N-1,0:M-1,abs(C_S))
            xlabel('time samples (Ts)')
            ylabel('Freq.samples')
            zlabel('C')
            title('spectrum method');
            subtitle(['mesh graph: L=' num2str(L) ', Kc=' int2str(Kc_value(i)) ', fdTs=' num2str(fdTs(j))])
            subplot(1,2,2)
            surf(1:N, 0:M-1,abs(C_S), 'MeshStyle', 'row'), view(0,90)
            xlabel('time samples (Ts)')
            ylabel('Freq.samples')
            subtitle(['surf graph: L=' num2str(L) ', Kc=' int2str(Kc_value(i)) ', fdTs=' num2str(fdTs(j))])
        end
    end
end

% Filter method
function c_f = filtermethod(fD,Ts,Ns,Kc)
N = Ns/10;
t = (-N*Ts:Ts:N*Ts);
g = besselj(1/4,2*pi*fD*abs(t))./nthroot(abs(t),4);
g(t==0) = nthroot(pi*fD,4)/gamma(5/4);
g=g/norm(g);
x = (randn(Ns,1) + 1i*randn(Ns,1))/sqrt(2);
c_f = conv(x,g) + Kc; 
end

% Spectrum method
function c_s = spectrummethod(fD,Ts,Ns,Kc)
fs = 1/Ts;
f = 0:fs/Ns:fD;
Sc = (1/(pi*fD)).*(1./sqrt(1-(f./fD).^2));
G = sqrt(Sc);
Gp = zeros(1,Ns);
Gp(1:length(G)) = G;
Gp(end-length(G)+1:end) = flip(G);
a = sqrt(Ns^2/norm(Gp)^2);
X = (randn(Ns,1) + 1i*randn(Ns,1))*a/sqrt(2); 
C = X'.*Gp;
c_s = ifft(C)+Kc;
end