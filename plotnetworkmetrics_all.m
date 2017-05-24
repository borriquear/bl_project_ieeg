function [netmetlist] = plotnetworkmetrics_all(Badjall, wiring_matrices, ie_patcondfreq, binary, typeofconnectivity)
%% plotnetworkmetrics plot network metrics 
%INPUT: Badjall, wiring_matrices, ie_patcondfreq, binary
% global nbofcalledfrequnecies;
% global totalfrequencies;
% totalfrequencies = ie_patcondfreq(end) - ie_patcondfreq(end-1) +1;
% nbofcalledfrequnecies = 1;
global cufreq;
global nbfreqs;
global lengend_msg;
lengend_msg = {'delta', 'theta', 'alpha', 'low beta', 'high beta', 'gamma'}; 
[nbpats, nbconds, nbfreqs, two]= size(Badjall);

if nbconds == 3
    %get only EC HYP
    cond1 = 'EC'; cond2 = 'HYP'
    condloop = 1:2:3; %Ec, HYP
else
    cond1 = 'EC'; cond2 = 'EO';
    condloop = 1:2; %EC EO
end
initpat = ie_patcondfreq(1); 
%initfreq = ie_patcondfreq(end); 
freqlist = wiring_matrices.frequencylist; 
freqlist = freqlist(ie_patcondfreq(end-1):ie_patcondfreq(end));

netmetlist = cell(nbpats, nbconds, nbfreqs);

%nbfreqs = length(freqlist);
listofpats = {}; listoffrequencies=[];
for pati=1:nbpats
    cupat = wiring_matrices.patientslist(pati + initpat - 1);
    %for plotting the patientid in the x-axis
    listofpats(pati) = cupat;
    %hdlf = figure;
    %hdlsolo = figure;
    for condj=condloop
        cucond =  wiring_matrices.conditionslist(condj);
        for freqk=1:nbfreqs
            cufreq = freqlist(freqk);
            listoffrequencies(freqk) = cufreq;
            netmetlist_el = {};
            %hdlf = figure;
            Badj = Badjall{pati, condj, freqk,1};
            thresholdv = Badjall{pati, condj, freqk,2};
            if condj == 2
                %plotting the contrast betwee two condition, plot - for EC and -- EO 
                linespec = '--';
            else
                linespec = '-';
            end
            netmetlist_el = plotnetworkmetrics(Badj, thresholdv, linespec, binary);
            netmetlist{pati, condj, freqk} = netmetlist_el;
            hold on
            %plotnetworkmetrics(handler, Badj, thresholdv)
        end
    end
    %figure(hdlf)
    if binary == 1
        msgtitle = sprintf('Network metrics Binarized Wiring Cost %s-based %s(-)-%s(--) pat=%s',typeofconnectivity,cond1, cond2, cupat{1})
    else
        msgtitle = sprintf('Wiring Cost %s-based %s (-) %s (--), pat=%s',typeofconnectivity, cond1, cond2,cupat{1})
    end
    title(msgtitle);
end

%calculate the difference between the two conditions
%mencontrast_fhdl = figure;
netmetlist_contrast_all= {};
typeofmetric = 3; % 3 wiring cost 1 for edges counter
for pati=1:nbpats
    for freqk=1:nbfreqs
        %eyes closed - eyes open
        cond1_v = netmetlist{pati, condloop(1), freqk}; 
        cond2_v= netmetlist{pati, condloop(end), freqk};
        netmetlist_contrast_all{pati,freqk} =cond1_v{typeofmetric} - cond2_v{typeofmetric} ;
        meanvectorofwc(pati,freqk) = mean(netmetlist_contrast_all{pati,freqk});
        stdvectorofwc(pati,freqk) = std(netmetlist_contrast_all{pati,freqk});
    end
end

%plot for mean all patients in frequency
figure;
maxnbofnetwors = 0; % the number of networks is variable, electrodes*electrodes
for i=1:nbpats
    plot(netmetlist_contrast_all{i,3})
    if length(netmetlist_contrast_all{i,3}) > maxnbofnetwors 
        maxnbofnetwors = length(netmetlist_contrast_all{i,3})
    end
    hold on
end
msgtitle = sprintf('mesoscopic Wiring Cost Difference(%s-%s) per subject for all networks %s-based in %s band', cond1, cond2, typeofconnectivity, 'alpha')
xlabel('network(threshold)'), ylabel('mean wiring cost difference each network');
xlim([0 maxnbofnetwors])
title(msgtitle)

% printeachpat = 0;
% if printeachpat == 1
%     plots_perfreq = {};
%     figure;
%     for freqk=1:nbfreqs
%         figure
%         bar(meanvectorofwc(:,freqk))
%         %errorbar(meanvectorofwc,stdvectorofwc)
%         msgtitle = sprintf('Mean Wiring Cost Difference(%s-%s) for all %s-based Freq=%.2f', cond1, cond2, typeofconnectivity, freqlist(freqk))
%         set(gca,'XTickLabel',listofpats(1:end));
%         xlabel('patients'), ylabel('mean wiring cost difference');
%         title(msgtitle)
%     end
% end
% %plot mean for each frequency, only 1 figure
% figure
% for i=1:length(listoffrequencies)
%     [freq_band] = getgreeksymbolfreq(listoffrequencies(i))
%     listoffrequencies_greek{i} = freq_band;
% end

%plot the mean differences for WC, for one patient and 1 frequency band(,alpha=3)
% figure
% plot(netmetlist_contrast_all{6,3})
% xlabel('network'), ylabel(' wiring cost difference ec-eo');
% title('Wiring cost difference EC-EO in Alpha band')
%plot the mean differences for number of edges, for one patient and 1 frequency band(,alpha=3)

end

function [metricsall] = plotnetworkmetrics(BM, delta, linespec, binary)
%% plotnetworkmetrics
%INPUT: hdl figure handler, Badj{1,:} Binary matrices, one for each threshold
%Badj{2,:} network metrics, one for each threshold, delta is the vector of
%threshold values, linespec = '--' eo  linespec = '-' ec
%OUTPUT: metricsall
%metricsall_contrast = {};
%global nbofcalledfrequnecies;
if binary == 0
    fprintf('plot wiringcost for weighted matrices\n')
    metricsall = plotwiringcostforweightedmatrices(hdl, BM, delta, linespec);
    %nbofcalledfrequnecies = nbofcalledfrequnecies +1;
    return
else
    [fi, col] = size(BM);
    netmetlist = keys(BM{2,1});
    for i=1:col
        %'B0'    'clustering'    'density'    'pathlength'
        mapobj = BM{2,i};
        Betti_v(i) = mapobj(netmetlist{1});
        clustering_v(i) = mapobj(netmetlist{2});
        wiringcost_v(i) = mapobj(netmetlist{3});
        pathlength_v(i) = mapobj(netmetlist{4});
    end
%     figure(hdl)
%     set(gca, 'XTick',[min(delta) mean(delta) max(delta)]);
%     %ylim([0 1]);
%     xlabel('delta'), ylabel('network metrics');
%     msgt = sprintf('Network metrics for BN threshold');
    %Change this to plot not only WC
    if 1==1 
            if strcmp(linespec, '--')
                plot(delta,wiringcost_v,'Color','b','LineStyle','-');hold on;
            else
                plot(delta,wiringcost_v,'Color','r','LineStyle','-');hold on;
            end
    else
%         plot(delta,Betti_v,'Color','m','LineStyle',linespec);hold on;
%         plot(delta,clustering_v,'Color','r','LineStyle',linespec);hold on;
%         plot(delta,pathlength_v,'Color','g','LineStyle',linespec);hold on;
%         % plot(delta,transitivity_v);
%         legend('Betti 0','clustering','wiringcost','pathlength')
    end
    metricsall = {Betti_v, clustering_v, wiringcost_v, pathlength_v};
end
end

function wiringcost_met = plotwiringcostforweightedmatrices(hdl, BM, delta, linespec)
%%plotwiringcostforweightedmtrices plot wiringcost y-axis, delta in x-axis
% global nbofcalledfrequnecies;
% global totalfrequencies;
global cufreq;
global nbfreqs;
global lengend_msg;
[fi, col] = size(BM);
netmetlist = keys(BM{2,1});
wiringcost =zeros(1,col);
for i=1:col
    mapobj = BM{2,i};
    wiringcost(i) = mapobj(netmetlist{1});
    %for other parameters put here
end
frequencycolors = {'y','m','c','r','g','b','w','k'};
%1.0000    1.7487    3.0579    5.3472    9.3506   16.3512   28.5930   50.0000
if cufreq >3 && cufreq < 4
    fqcolor = 1; %yellow for delta
elseif cufreq > 5 && cufreq < 6
    fqcolor = 2; %m for theta
elseif cufreq > 8 && cufreq < 10
    fqcolor = 3; %c for alpha
elseif cufreq > 10 && cufreq < 20
    fqcolor = 4; %r for low beta
elseif cufreq > 20 && cufreq < 30
    fqcolor = 5; %g for high beta
elseif cufreq > 30 && cufreq < 55
    fqcolor = 6; %b for gamma
end
%fqcolor = rem(nbofcalledfrequnecies,totalfrequencies) +div(nbofcalledfrequnecies,totalfrequencies);
%y yello m magenta c cyan r red g green b blue w white k black

% figure(hdl)
% set(gca, 'XTick',[min(delta) mean(delta) max(delta)])
% %ylim([0 1]);
% xlabel('threshold'), ylabel('wiring cost')
% %msgt = sprintf('Wiring cost EC(-) EO(--) weighted networks \');
% %legend(conditionslist{1}, conditionslist{2}, 'Interpreter', 'None')
% %title(msgt,'interpreter', 'none');
% 
% %plot(delta,wiringcost,'Color','m','LineStyle',linespec);hold on;
% plot(delta,wiringcost,'Color',frequencycolors{fqcolor},'LineStyle',linespec);hold on;
% legend(lengend_msg(1:nbfreqs) )
% wiringcost_met = {wiringcost};
%%%
end