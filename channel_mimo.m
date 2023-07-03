function chan = channel_mimo(cmode,L,rx,tx,severity)
% This is a mimo channel generator for M-Malaga and double gamma
% channels, based on:  
% "A. Jurado-Navas, J. Maria, J. Francisco, and A. Puerta-Notario, “A unifying statistical
% model for atmospheric optical scintillation,” in Numerical Simulations of Physical and
% Engineering Processes. InTech, sep 2011 
% and 
% K. Yiannopoulos, N. C. Sagias, and A. C. Boucouvalas, “Average error probability of an
% optically pre-amplified pulse-position modulation multichannel receiver under malaga-m
% fading,” Applied Sciences, vol. 10, no. 3, 2020 
%
% Created for the A.Aspreas and K. Yiannopoulos, Reed Solomon Error Correction in Pre-amplified Pulse
% Position Modulation Receivers"  AEUE - International Journal of Electronics and Communications 
%
% Created at 03/07/2023 by Antonios Aspreas

%% Inputs:
% cmode: Channel type. 
%   1. No Channel
%   2. Gamms-gamma channel 
%   3. M-Malaga channel 
% L: Number of channel amplitudes 
% rx: Number of receivers
% tx: Number of transmitters
% severity: Weak, Medium,Strong
%
% Example chan = channel_mimo(3,1e6,1,1,'weak');

if ~isnumeric(cmode) || ~rem(cmode,1)==0
    error('Channel type must be an integer between 1-3');
end
if ~isnumeric(L)|| ~rem(L,1)==0
    error('Number of samples L must be an integer');
end
if ~isnumeric(rx)|| ~rem(rx,1)==0
    error('Number of receivers rx must be an integer');
end
if ~isnumeric(tx)|| ~rem(tx,1)==0
    error('Number of transmitters tx must be an integer');
end
if ~ischar(severity) && ~isstring(severity)
    error('Severity type must be a string');
end
switch cmode
    case 1
        chan = 1;
    case 2
        switch severity
            case 'weak'
                %1-100m
                a = 16.5347;
                b = 14.9057;
            case 'strong'
                % 1-350m
                a = 4.22772;
                b = 2.3177;
            otherwise
                error('No case %s for gamma-gamma channel. Choose between weak,strong',severity)
        end
        chan = gamrnd(a,1,L,rx,tx).*gamrnd(b,1,L,rx,tx)/(a*b*rx*tx);
    case 3
        switch severity
            case 'weak'
                a=50;
                b = 14;
                omega=1.0621;
                beta0 = 0.0216;
                rho = 0.86; 
                deltaphi=0;
            case'moderate'
                a=2.55;
                b = 22;
                omega=0.4618;
                beta0 = 0.6525;
                rho = 0.988; 
                deltaphi=pi/2;
            case'strong'
                a=2.2814;
                b = 33;
                omega=1.33;
                beta0 = 0.4231;
                rho = 0.84; 
                deltaphi=0;
            otherwise
                error('No case %s for m-malaga channel. Choose between weak,moderate,strong',severity)
        end
        gam = 2*beta0*(1-rho);
        omegaprime = omega+rho*2*beta0 +2*sqrt(2*beta0*rho)*cos(deltaphi);

        g = sqrt(gamrnd(b,1,L,rx,tx)/b);
        USprime = sqrt(2*beta0)*randn(L,rx,tx);
        phi = rand(L,rx,tx)*2*pi;
        r = g.*(sqrt(omega)+sqrt(2*beta0*rho).*exp(1i*deltaphi) +sqrt(1-rho).*USprime.*exp(1i*phi) );
        
        y = abs(r.^2);
        x = gamrnd(a,1,L,rx,tx)/a;
        irradiance = x.*y;
        f_mean=malaga_moment(1,a,b,gam,omegaprime);
        %         f_std=malaga_moment(2,a,b,gam,omegaprime);
        
        chan = (irradiance /f_mean );
        %         chan = irradiance ;
        
    otherwise
        error('No channel %s. \n Choose between 1: No channel, 2:gamma-gamma, 3: M-Malaga',num2str(cmode))
end


end