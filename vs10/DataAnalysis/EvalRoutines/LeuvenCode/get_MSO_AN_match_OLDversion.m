clc;clear;

MSOdatapath = 'C:\Users\Mark\Dropbox\BinRev';
ANdatapath = 'C:\Users\Mark\Dropbox\MonRev';

ANdatalist = load('C:\Users\Mark\Dropbox\Disparity_MATLAB\MonRevData.mat');

MSOfolderlist = dir(fullfile(MSOdatapath,'*_*'));
nrMSOunits = length(MSOfolderlist);

% for i = 1:nrMSOunits
for i = 1:nrMSOunits
    MSOunitID = MSOfolderlist(i).name;
    if ~isempty(dir(fullfile(MSOdatapath,MSOunitID,'*BinRev_70.00dB.mat')))
        load(fullfile(MSOdatapath,MSOunitID,[MSOunitID '_BinRev_70.00dB.mat']));
        for j = 1:2
            DF = kernels.h1cdomfreq(j);
            lowF = DF*2^-0.075;
            highF = DF*2^0.075;
            unitinds = find(ANdatalist.DF>=lowF & ANdatalist.DF<=highF & ANdatalist.CF<=2500);
            nrpossunits = length(unitinds);
            MSOphase = kernels.h1cphase(:,j);
            MSOweights = kernels.h1coeffs.Monaural.weights(:,j);
            Freq = kernels.h1ffax;
            %Get all time and phase delays for this ear between AN and MSO
            for k = 1:nrpossunits
                ANunitID{j,k} = ANdatalist.unitID{unitinds(k)};
                AN = load(fullfile(ANdatapath,ANunitID{j,k},[ANunitID{j,k} '_MonRev_70.00dB.mat']));
                ANphase{j}(:,k) = AN.kernels.h1cphase;
                ANweights{j}(:,k) = AN.kernels.h1coeffs.Monaural.weights;
                W{j}(:,k) = ANweights{j}(:,k).*MSOweights;
                W{j}(:,k) = W{j}(:,k)./max(W{j}(:,k));
                AN_MSO_Phase{j}(:,k) = (MSOphase-ANphase{j}(:,k))/(2*pi);
                [AN_MSO_PhaseDelay{j}(k),AN_MSO_TimeDelay{j}(k)] = fit_circ_lin(Freq(W{j}(:,k)>0),AN_MSO_Phase{j}(W{j}(:,k)>0,k),W{j}(W{j}(:,k)>0,k));
            end
        end
        %Compare AN phase contra vs ipsi to get the cochlear CD and
        %cochlear CP for this unit
        nrANunitsIPSI = size(ANphase{1},2);
        nrANunitsCONTRA = size(ANphase{2},2);
        [ipsiANinds,contraANinds] = meshgrid(1:nrANunitsIPSI,1:nrANunitsCONTRA);
        ipsiANinds = ipsiANinds(:);
        contraANinds = contraANinds(:);
        nrcomps = length(contraANinds);%total number of comparisons possible for these pairs
        for j = 1:nrcomps
            XANW(:,j) = ANweights{1}(:,ipsiANinds(j)).*ANweights{2}(:,contraANinds(j));%Cross nerve weights
            XANW(:,j) = XANW(:,j)/max(XANW(:,j));
            XANP(:,j) = (ANphase{2}(:,contraANinds(j))-ANphase{1}(:,ipsiANinds(j)))/(2*pi);%Cross nerve phase
            [XAN_PhaseDelay(j),XAN_TimeDelay(j)] = fit_circ_lin(Freq(XANW(:,j)>0),XANP(XANW(:,j)>0,j),XANW(XANW(:,j)>0,j));
        end
        CD_cochlea.mean = mean(XAN_TimeDelay);
        CD_cochlea.std = std(XAN_TimeDelay);
        CD_cochlea.sem = CD_cochlea.std/(sqrt(nrcomps));
        
        [CP_cochlea.mean,R] = circmean(XAN_PhaseDelay*(2*pi)); %returns CP_cochlea.mean in radians on the interval [0, 2pi]
        CP_cochlea.mean = wrapToPi(CP_cochlea.mean)/(2*pi); %Put it back into cycles on the interval [-0.5, 0.5]
        CP_cochlea.std = sqrt(log(1/R^2)); %circular standard deviation
        CP_cochlea.sem = CP_cochlea.std/sqrt(nrcomps);
        
        
        
        %Compare the AN to MSO phase contra vs ipsi to get the
        %brainstem CD and brainstem CP.
        nrunitsIPSI = size(AN_MSO_Phase{1},2);
        nrunitsCONTRA = size(AN_MSO_Phase{2},2);
        [ipsiinds,contrainds] = meshgrid(1:nrunitsIPSI,1:nrunitsCONTRA);
        ipsiinds = ipsiinds(:);
        contrainds = contrainds(:);
        nrcomps = length(contrainds);%total number of comparisons possible for these pairs
        for j = 1:nrcomps
            XBW(:,j) = W{1}(:,ipsiinds(j)).*W{2}(:,contrainds(j));%Cross brain weights
            XBW(:,j) = XBW(:,j)/max(XBW(:,j));
            XBP(:,j) = AN_MSO_Phase{2}(:,contrainds(j))-AN_MSO_Phase{1}(:,ipsiinds(j));%Cross brain phase
            [XB_PhaseDelay(j),XB_TimeDelay(j)] = fit_circ_lin(Freq(XBW(:,j)>0),XBP(XBW(:,j)>0,j),XBW(XBW(:,j)>0,j));
        end
        CD_brain.mean = mean(XB_TimeDelay);
        CD_brain.std = std(XB_TimeDelay);
        CD_brain.sem = CD_brain.std/(sqrt(nrcomps));
        
        [CP_brain.mean,R] = circmean(XB_PhaseDelay*(2*pi)); %returns CP_brain.mean in radians on the interval [0, 2pi]
        CP_brain.mean = wrapToPi(CP_brain.mean)/(2*pi); %Put it back into cycles on the interval [-0.5, 0.5]
        CP_brain.std = sqrt(log(1/R^2)); %circular standard deviation
        CP_brain.sem = CP_brain.std/sqrt(nrcomps);
        
        
        save(fullfile(MSOdatapath,MSOunitID,[MSOunitID '_AN_MSO_match_70.00dB.mat']),'ANunitID','ANphase','ANweights',...
            'W','AN_MSO_Phase','AN_MSO_PhaseDelay','AN_MSO_TimeDelay','XANW','XANP','XAN_PhaseDelay',...
            'XAN_TimeDelay','CD_cochlea','CP_cochlea','XBW','XBP','XB_PhaseDelay','XB_TimeDelay','CD_brain','CP_brain');
        
        clear ('ANunitID','ANphase','ANweights',...
            'W','AN_MSO_Phase','AN_MSO_PhaseDelay',...
            'AN_MSO_TimeDelay','XANW','XANP',...
            'XAN_PhaseDelay','XAN_TimeDelay','CD_cochlea',...
            'CP_cochlea','XBW','XBP','XB_PhaseDelay',...
            'XB_TimeDelay','CD_brain','CP_brain');
    end
end
