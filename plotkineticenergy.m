function kinetic_matrices = plotkineticenergy(wiring_matrices, powerspecmatrix_freqbands)
%plotkineticenergy points per patient (x,y) x is the percentage of power relative 
%to the other bands and y is the K for that band
%plot bars of kinetic energy (P*F*frequency)^2

% 1/2v^2 is the kinetic energy to transport a punctual mass m=1 between two
% points, electrode positions.
%IN wiring_matrices, powerspecmatrix_freqbands

global nbpats;
global freqlist;
global nbfreqs;
global nbconds;
globalFsDir = loadglobalFsDir();
kinetic_matrices = struct ;
nbpats = size(wiring_matrices.patientslist,2);
freqlist = wiring_matrices.frequencylist;
nbfreqs = length(freqlist);
nbconds = size(wiring_matrices.conditionslist, 2); nbconds = 2;

%nbfreqs = length( wiring_matrices.frequencylist);
barmatrix_ispc = zeros(nbpats, nbconds*nbfreqs);
barmatrix_pli = zeros(nbpats, nbconds*nbfreqs);
barmatrix_power = zeros(nbpats, nbconds*nbfreqs);
barmatrix_chg_ispc = zeros(nbpats, nbconds*nbfreqs);
barmatrix_chg_pli = zeros(nbpats, nbconds*nbfreqs);
barmatrix_chg_power = zeros(nbpats, nbconds*nbfreqs);
%calculate the kinetic matrices
initfreq= 1;
wiring_ispc ={};
%wiring_ispc = p*f
%fm factor from mm to m
fm = 10^-3;
for p=1:nbpats
    for ic=1:nbconds
        for f=initfreq:nbfreqs
            wiring_ispc{p,ic,f} = wiring_matrices.wiring_ispc{p,ic,f}*wiring_matrices.frequencylist(f);
            wiring_ispc{p,ic,f}= 0.5*(wiring_ispc{p,ic,f}*fm)^2;
            wiring_pli{p,ic,f} = wiring_matrices.wiring_pli{p,ic,f}*wiring_matrices.frequencylist(f);
            wiring_pli{p,ic,f}= 0.5*(wiring_pli{p,ic,f}*fm)^2;
            wiring_power{p,ic,f} = wiring_matrices.wiring_power{p,ic,f}*wiring_matrices.frequencylist(f);
            wiring_power{p,ic,f}= 0.5*(wiring_power{p,ic,f}*fm)^2;
            wiring_distperfreq{p,ic,f}= 0.5*(wiring_matrices.distMatrixcell{p}*f)^2;
            
        end
    end
end
kinetic_matrices.wiring_ispc = wiring_ispc;
kinetic_matrices.wiring_pli = wiring_pli;
kinetic_matrices.wiring_power = wiring_power;
kinetic_matrices.patientslist = wiring_matrices.patientslist;
kinetic_matrices.conditionslist = wiring_matrices.conditionslist;
kinetic_matrices.frequencylist = wiring_matrices.frequencylist;
kinetic_matrices.distMatrixcell = wiring_matrices.distMatrixcell;
kinetic_matrices.wiring_distperfreq = wiring_distperfreq;
%save the kinetic_matrices structure
matfilename = 'wiringcost_matrices.mat';
matfilename = fullfile(globalFsDir, matfilename);
disp(kinetic_matrices)
save(matfilename,'kinetic_matrices');

%build the power spectrum matrix PS
initfreq = 3;
PS = zeros(nbconds,nbpats,nbfreqs- initfreq);
for ic=1:nbconds
    for inp=1:nbpats
        ele = powerspecmatrix_freqbands(inp,ic); ele = ele{1};
        for ifq=1:nbfreqs-initfreq
            PS(ic, inp, ifq ) = mean(ele(:,ifq));
        end
    end
end
%Build the Kietic Energy matrix for the different funtional connectivity
%measures
K_ispc = zeros(nbconds,nbpats,nbfreqs - initfreq); 
K_pli = zeros(nbconds,nbpats,nbfreqs - initfreq);
K_power = zeros(nbconds,nbpats,nbfreqs - initfreq);
for ic=1:nbconds
    for inp=1:nbpats
        for ifq=1:nbfreqs-initfreq
            ele= kinetic_matrices.wiring_ispc(inp,ic,ifq); 
            K_ispc(ic,inp, ifq) = mean(ele{1}(:));
            ele= kinetic_matrices.wiring_pli(inp,ic,ifq); 
            K_pli(ic,inp, ifq) = mean(ele{1}(:));
%             ele= kinetic_matrices.wiring_power(inp,ic,ifq); 
%             K_power(ic,inp, ifq) = mean(ele{1}(:));
ele = kinetic_matrices.wiring_distperfreq(inp,ic,ifq); 
 K_power(ic,inp, ifq) = mean(ele{1}(:));

        end
    end
end

%plot two vectors, x power per band
%y kinetic energy, both are normalized
condlis = {'EC' , 'EO '};
for ic=1:nbconds
    %one figure per condition, nbfreqs data points per patient
    figure
    y_ispcT = []; y_pliT = []; y_powerT = [];xT =[];
    for ip=1:nbpats
        x = zeros(1, nbfreqs-initfreq); 
        y_ispc = zeros(1, nbfreqs-initfreq);
        y_pli= zeros(1, nbfreqs-initfreq);
        y_power= zeros(1, nbfreqs-initfreq);
        for ifq=1:nbfreqs-initfreq
            y_ispc(ifq) = K_ispc(ic,ip, ifq);
            y_pli(ifq) = K_pli(ic,ip, ifq);
            y_power(ifq) = K_power(ic,ip, ifq);
            x(ifq) = PS(ic,ip, ifq);
        end
        %plot x,y nbfreqs-initfreq data points
        sumx= sum(x);
        sumy_ispc = sum(y_ispc);
        x = x/sumx; y_ispc = y_ispc/sumy_ispc;
        xT = [xT x];
        R_ispc= corrcoef(x,y_ispc); R_ispc_sq = R_ispc(2)^2;
        vR_ispc(ip) = R_ispc(2); vR_ispc_sq(ip) = R_ispc_sq;
        subplot(1,3,1)
        y_ispcT= [y_ispcT y_ispc]
        subplot(1,3,2)
        sumy_pli = sum(y_pli);
        y_pli = y_pli/sumy_pli;
        R_pli = corrcoef(x,y_pli); R_pli_sq = R_pli(2)^2;
        vR_pli(ip) = R_pli(2); vR_pli_sq(ip) = R_pli_sq;
        y_pliT = [y_pliT y_pli]
        subplot(1,3,3)
        sumy_power = sum(y_power);
        y_power = y_power/sumy_power;
        R_power = corrcoef(x,y_power); R_power_sq = R_power(2)^2;
        vR_power(ip) = R_power(2); vR_power_sq(ip) = R_power_sq;
        y_powerT = [y_powerT y_power]
    end  
    disp([vR_ispc vR_pli vR_power(ip)])
    
    subplot(1,3,1)
    %Fit with a two term exponential model
    f = fit(xT',y_ispcT', 'exp2')
    plot(f,xT,y_ispcT)
    xlabel('mean power per band'),ylabel('K= 1/2 (P*w_R*fq)^2')
    text(0.3,max(y_ispcT)/2, ['R=', num2str(mean(vR_ispc))]);
    
    subplot(1,3,2)
    f = fit(xT',y_pliT', 'exp2')
    plot(f,xT,y_pliT)
    xlabel('mean power per band'),ylabel('K= 1/2 (P*w_{PLI}*fq)^2')
    text(0.3,max(y_pliT)/2, ['R=', num2str(mean(vR_pli))]);
    stitle = sprintf('Kinetic Energy, %s',condlis{ic})
    title(stitle);
    
%     subplot(1,3,3)
%     f = fit(xT',y_powerT', 'exp2')
%     plot(f,xT,y_powerT)
%     xlabel('mean power per band'),ylabel('K= 1/2 (P*WPw*fq)^2')
%     text(0.3,max(y_powerT)/2, ['R=', num2str(mean(vR_power))]);
    
    subplot(1,3,3)
    f = fit(xT',y_powerT', 'exp2')
    plot(f,xT,y_powerT)
    xlabel('mean power per band'),ylabel('K= 1/2 (P*WPw*fq)^2')
    text(0.3,max(y_powerT)/2, ['R=', num2str(mean(vR_power))]);
    
end
end