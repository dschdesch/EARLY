function varargout = getLewisSTRF(varargin)
%Returns STRF from 2nd-order Wiener kernel.
%Implements the algorithm described by Lewis & van Dijk (2004), as modified
%by Sneary & Lewis (2007).
%--------------------------------------------------------------------------
%Mark Sayles. 2015. Leuven.
%--------------------------------------------------------------------------
switch nargin
    case 2
        dataIN = varargin{1};
        SR = varargin{2};
        flag = 0;
    case 3
        dataIN = varargin{1};
        SR = varargin{2};
        NulldataIN = varargin{3};
        flag = 0;
    case 4
        dataIN = varargin{1};
        SR = varargin{2};
        NulldataIN = varargin{3};
        flag = varargin{4};
end
nT = size(dataIN,1);
nFFT = 2^nextpow2(2*nT-1);
dataINflip = fliplr(dataIN);
S = zeros(nFFT,nT);
LewisFreq = (SR/4*linspace(0,1,nFFT/4));
LewisTime = (1/SR)*(0:nT-1);
for T = 1:nT
    d1 = diag(dataINflip,nT-(2*T)+1);
    d2 = diag(dataINflip,nT-(2*T));
    nd = length([d1;d2]);
    d = zeros(1,nd);
    if (2*T)<=nT
        d(1:2:end) = d2;
        d(2:2:end) = d1;
    else
        d(1:2:end) = d1;
        d(2:2:end) = d2;
    end
    d = [d(ceil(nd/2):end) d(1:floor(nd/2))];
    S(1:ceil(nd/2),T) = d(1:ceil(nd/2));
    S(nFFT-floor(nd/2)+1:nFFT,T) = d(ceil(nd/2)+1:nd);
end
G = real(fft(S));
G = G(1:nFFT/4,:);%STRF to SR/4
% G(LewisFreq<200,:) = 0;
LSTRF = G;
% SmoothSTRF = smooth2a(G,5,5);
%Get the modulation tuning function and modulation transfer function
% [ModFn,MTF,ModScale.y] = get_STRF_MTF(LSTRF,LewisFreq(2),nFFT);
% ModScale.x = (SR/4)*linspace(0,1,nT/4+1);
% ModScale.x = ModScale.x(2:end);%No 0-Hz component

%--------------------------------------------------------------------------
%If we have the null data define some statistically significant portions of the STRF
%--------------------------------------------------------------------------
if exist('NulldataIN','var')
    for i=1:size(NulldataIN,3)
        NulldataINflip = fliplr(NulldataIN(:,:,i));
        S = zeros(nFFT,nT);
        for T = 1:nT
            d1 = diag(NulldataINflip,nT-(2*T)+1);
            d2 = diag(NulldataINflip,nT-(2*T));
            nd = length([d1;d2]);
            d = zeros(1,nd);
            if (2*T)<=nT
                d(1:2:end) = d2;
                d(2:2:end) = d1;
            else
                d(1:2:end) = d1;
                d(2:2:end) = d2;
            end
            d = [d(ceil(nd/2):end) d(1:floor(nd/2))];
            S(1:ceil(nd/2),T) = d(1:ceil(nd/2));
            S(nFFT-floor(nd/2)+1:nFFT,T) = d(ceil(nd/2)+1:nd);
        end
        temp = real(fft(S));
%         temp(LewisFreq<200,:) = 0;
        LSTRFNull(:,:,i) = temp(1:nFFT/4,:);%STRF to SR/4;
        clear temp;
    end
%     STRFstd = std(LSTRFNull(:));
    ZSTRF = LSTRF./std(LSTRFNull,0,3);
    %Find all "significant" regions
    z_crit1 = 3;
    z_crit2 = 9;
%     areacrit = 0.1;
    areacrit = 1; %Minimum area criterion
    z_crit1 = [-1*z_crit1 z_crit1]; %Z-score criteria
    figure(999);
    [dummy,h] = contour3(LewisTime*1000,LewisFreq/1000,smooth2a(ZSTRF,5,5),z_crit1);
    ex_count = 0;
    in_count = 0;
    for j=1:length(h)
        zz=get(h(j),'Zdata');
        x=get(h(j),'XData');
        y=get(h(j),'YData');
        x(end)=x(1);
        y(end)=y(1);
        ar=polyarea(x,y);
        Ind = inpolygon(repmat(LewisTime*1000,numel(LewisFreq),1),repmat(LewisFreq'/1000,1,numel(LewisTime)),x,y);
        if zz(1)<0
            if ar>areacrit || (~isempty(max(max(abs(ZSTRF(Ind))))) && max(max(abs(ZSTRF(Ind))))>=z_crit2);
                in_count = in_count+1;
                suppression_curve{in_count}.x = x;
                suppression_curve{in_count}.y = y;
            end
        elseif zz(1)>0
            if ar>areacrit || (~isempty(max(max(abs(ZSTRF(Ind))))) && max(max(abs(ZSTRF(Ind))))>=z_crit2);
                ex_count = ex_count+1;
                excitation_curve{ex_count}.x = x;
                excitation_curve{ex_count}.y = y;
            end
        end
    end
end
close(999);
if ~exist('excitation_curve','var')
    excitation_curve = [];
end
if ~exist('suppression_curve','var')
    suppression_curve = [];
end

%--------------------------------------------------------------------------
%Return the output
%--------------------------------------------------------------------------

varargout{1} = LSTRF;
varargout{2} = LewisFreq;
varargout{3} = LewisTime;
if exist('NulldataIN','var')
    varargout{4} = ZSTRF;
    varargout{5} = excitation_curve;
    varargout{6} = suppression_curve;
end
% varargout{1} = LSTRF;
% varargout{2} = LewisFreq;
% varargout{3} = LewisTime;
% varargout{4} = ModFn;
% varargout{5} = MTF;
% varargout{6} = ModScale;
% if exist('NulldataIN','var')
%     varargout{7} = ZSTRF;
%     varargout{8} = excitation_curve;
%     varargout{9} = suppression_curve;
%     varargout{10} = MTFNull;
% end

%Locals
    function [filtout] = local_high_pass(IN,cut1,cut2)
        fcuts = [cut1 cut2];
        fsamp = SR; %sample frequency
        mags = [0 1];
        devs = [0.01 10^(-40/20)];% 1% passband ripple and 40-dB attenuation of stopband
        [n,Wn,beta,ftype] = kaiserord(fcuts,mags,devs,fsamp);
        b = fir1(n,Wn,ftype,kaiser(n+1,beta),'noscale'); %Filter coefficients to match the required parameters
        filtout = filtfilt(b,1,IN); %Zero-phase filtering
    end
end
