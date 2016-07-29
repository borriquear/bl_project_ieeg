function [ ] = savefigaspdfnomargins()
%savefigaspdfnomargins save a .fig as a pdf without huge margins
% Open only one figure at a time and change the pdfname prior to run this code

pdfname = 'power-amount-42-allcond-allbands.pdf';
path_fig_files = 'C:\workspace\github\figures';
filename = fullfile(path_fig_files,pdfname);
fh=findobj(0,'type','figure');
save2pdf(filename,fh);

% nbofopenfigs = length(fh);
% h = figure(nbofopenfigs);
% filename = fullfile(path_fig_files,pdfname);
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3) ;
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig,filename,'-dpdf')
% disp(filename)
end

