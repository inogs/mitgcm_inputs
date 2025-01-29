%%%% create bottom fluxes (no point sources) %%%%%
% creato a partire da make_bottom_fluxes_AZAL_HR_V4.m (Gaeta)
% confrontato con make_bottom_fluxes_AZAL_V5.m (ultima versione per il Lazio)

opengl neverselect
close all
clear all
curdir=pwd;

nx=450;
ny=300;
nz=37;
nt=365;

% Global Area (Integrated horizontal Area [m^2])
% GA=1.179E+09

delZ=[0.5;  0.5;  0.5;  0.5;  0.5;  0.5;  1.0;  1.0;  1.0;  1.0; ...
      1.0;  1.0;  1.0;  1.0;  1.0;  1.0;  1.0;  1.0;  1.0;  1.0; ...
      1.0;  1.0;  1.0;  1.0;  1.0;  1.0;  1.0;  1.0;  1.0;  1.0; ...
      1.0;  1.0;  1.0;  1.0;  1.0;  2.0;  2.0];

% Mean (±standard deviation) annual benthic fluxes of N, P, Si and O2 (mmol m−2 d−1).
% Reference                 NO−3		NH+4		PO3−4		Si(OH)4		O2
% Bertuzzi et al. (1997)	0.17±0.73	0.8±0.7		0.029±0.05	2.59±2.3	−20.4±8.9
% BFM–POM 1D				0.27±0.16	0.63±0.38	0.048±0.03	0.47±0.41	−5.14±3.5

Nb=(0.17+0.8)./86400; % [mmol/m^2/s]
Pb=0.029./86400; % [mmol/m^2/s]
Ob=-20.4./86400; % [mmol/m^2/s]

namef   =['N';'P';'B';'C';'O'];
patom   =[ 14; 31;  1;  1; 16];
bentflux=[ Nb; Pb; 0.; 0.; Ob];

load('mask_bottom_GoT_iNEST_V2'); % bottom
load('mask_qbottom_GoT_iNEST_V2'); % qbottom
load('mask_surface_GoT_iNEST_V2'); % surf

cd /share/data/squerin/esperimenti/GoT_iNEST/bathymetry

fidx=fopen('./bathydata_GoT_iNEST_V2/DXG.data','r','l');
dxg=fread(fidx,[nx,ny],'float32');
fclose(fidx);

fidy=fopen('./bathydata_GoT_iNEST_V2/DYG.data','r','l');
dyg=fread(fidy,[nx,ny],'float32');
fclose(fidy);

fidb=fopen('./bathydata_GoT_iNEST_V2/Depth.data','r','l');
btm=fread(fidb,[nx,ny],'float32');
fclose(fidb);

cd(curdir)

%%% FLUSSI %%%

% minlev=6; % 6(AZAL)=~15 m; 3.0010 (profondità minima)
% minlev=8; % 8(AZAL)=~25 m;
% minlev=28; % 28(AZAL_HR)=~25 m; 3.00 (profondità minima)
minlev=18; % 18(GoT_iNEST)=~15 m; 3.00 (profondità minima)
mindepth=sum(delZ(1:minlev));
% levlim=15; % -> ~70 m; levlim=12(AZAL) -> W-Interf. bottom 48.8900 (limite dei ~50 m)
% levlim=45; % -> ~70 m; levlim=42(AZAL_HR) -> limite dei ~50 m
levlim=37; % -> ~35 m; levlim=37(GoT_iNEST) -> limite dei ~35 m
exp=4;
mult=2;
lim=sum(delZ(1:levlim));
btml=btm;
btml(btml>=lim)=lim; % taglio il fondo a ~35 m
sp=mult.*((surf.*lim-btml)./(lim-mindepth)).^exp;
figure;pcolor(sp');shading flat;colorbar
spm=sp;
title(['minlev=' int2str(minlev) ' levlim=' int2str(levlim) ' exp=' int2str(exp)])
% caxis([0 1])

btml=bottom;
btml(btml>=levlim)=levlim; % taglio il fondo al livello levlim
sp2=mult.*((surf.*levlim-btml)./(levlim-minlev)).^exp;
figure;pcolor(sp2');shading flat;colorbar
title(['minlev=' int2str(minlev) ' levlim=' int2str(levlim) ' exp=' int2str(exp)])
% caxis([0 1])

for var=1:size(namef,1)
    eval(['fidout_',namef(var),...
        '=fopen(''',namef(var),'_bottom_fluxes_GoT_iNEST_V2_x2.dat'',''w'',''l'');']);
    tflux=zeros(1,nt);
    for t=1:nt
        cy(t)=0.5.*(1-cos(((2*pi./365).*t)+2*pi./365*(10+0))); %#ok<SAGROW>
        c1(t)=cy(t);
        % if t>=210
        %     cy(t)=cy(t).*1/(t-209).^0.5;
        %     c1(t)=c1(t);
        % end
        bflux=spm.*cy(t).*bentflux(var);
        eval(['fwrite(fidout_',namef(var),',bflux,''float32'');']);

        bf=bflux.*dxg.*dyg;
        bf(bf==0)=NaN;
        if any(~isnan(bf(:)))
            tflux(t)=nansum(nansum(bf)).*86400./1e+9.*patom(var); % [tonn/gg]
        else
            disp(tflux(t))
        end
    end
    yrflux=sum(tflux);
    disp(yrflux);
    figure;pcolor(bflux');shading flat;colorbar
    eval(['fclose(fidout_',namef(var),');'])
end
