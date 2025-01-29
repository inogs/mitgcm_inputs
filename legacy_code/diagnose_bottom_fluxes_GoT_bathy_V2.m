%%%% diagnose bottom fluxes (no point sources) %%%%%
opengl neverselect
close all
clear all
curdir=pwd;

nx=450;
ny=300;
% nz=50;
nt=365;

% Global Area (Integrated horizontal Area [m^2])
% GA=1.179E+09

% Mean (±standard deviation) annual benthic fluxes of N, P, Si and O2 (mmol m−2 d−1).
% Reference                 NO−3		NH+4		PO3−4		Si(OH)4		O2
% Bertuzzi et al. (1997)	0.17±0.73	0.8±0.7		0.029±0.05	2.59±2.3	−20.4±8.9
% BFM–POM 1D				0.27±0.16	0.63±0.38	0.048±0.03	0.47±0.41	−5.14±3.5

namef   =['N';'P';'B';'C';'O'];
patom   =[ 14; 31;  1;  1; 16];

load('mask_surface_GoT_iNEST_V2'); % surf
surf(surf==0)=NaN;

cd /share/data/squerin/esperimenti/GoT_iNEST/bathymetry

fidx=fopen('./bathydata_GoT_iNEST_V2/DXG.data','r','l');
dxg=fread(fidx,[nx,ny],'float32');
fclose(fidx);

fidy=fopen('./bathydata_GoT_iNEST_V2/DYG.data','r','l');
dyg=fread(fidy,[nx,ny],'float32');
fclose(fidy);

cd(curdir)

%%% FLUSSI %%%

for var=1:size(namef,1)
    if strcmp(namef(var),'N')
        eval(['fidin_',namef(var),...
            '=fopen(''',namef(var),'_bottom_fluxes_GoT_iNEST_V2_x2.dat'',''r'',''l'');']);
    elseif strcmp(namef(var),'P')
        eval(['fidin_',namef(var),...
            '=fopen(''',namef(var),'_bottom_fluxes_GoT_iNEST_V2_x3.dat'',''r'',''l'');']);
    else
        eval(['fidin_',namef(var),...
            '=fopen(''',namef(var),'_bottom_fluxes_GoT_iNEST_V2_flat.dat'',''r'',''l'');']);
    end
    
    tbflux=zeros(nx,ny,nt);
    tbf=zeros(nx,ny,nt);
    for t=1:nt
        eval(['bflux=fread(fidin_',namef(var),',[nx,ny],''float32'');']); % [mmol/m^2/s]
        tbflux(:,:,t)=bflux; % [mmol/m^2/s]
        tbf(:,:,t)=bflux.*dxg.*dyg; % [mmol/s]
    end
    tbfint=squeeze(sum(sum(tbf,1),2));
    
    tonnyr=sum(tbfint).*86400./1e+9.*patom(var); % [tonn/yr]
    disp([namef(var),': ',num2str(tonnyr,'%e'),' [tonn/yr]']);
    
    figure;plot(tbfint); % [mmol/s]
    title(['instant flux of ',namef(var),' - spatially integrated [mmol/s]'])
    set(gcf,'renderer','zbuffer','PaperPositionMode','auto')
    eval(['print -dpng -r600 instant_flux_',namef(var),'.png'])
    
    mbflux=mean(tbflux,3);
    figure;pcolor((mbflux.*surf)');shading flat;colorbar
    title(['mean flux of ',namef(var),' [mmol/m^2/s]'])
    set(gcf,'renderer','zbuffer','PaperPositionMode','auto')
    eval(['print -dpng -r600 mean_flux_',namef(var),'.png'])
    % da fare eventualmente plot mensili...
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure;
    for stat=20:20:240
        plot(squeeze(tbflux(200,stat,:)),'b')
        hold on
        plot(repmat(squeeze(mean(tbflux(165,stat,:),3)),365),'--r')
    end
    title(['mean bottom flux of ',namef(var),' = ',num2str(mean(mean(mbflux))),' [mmol/m^2/s]'])
    set(gcf,'renderer','zbuffer','PaperPositionMode','auto')
    eval(['print -dpng -r600 time_series_flux_transect_',namef(var),'.png'])
end
