%%%% create surface fluxes (atmospheric deposition) %%%%%
opengl neverselect
close all
clear all
curdir=pwd;

nx=450;
ny=300;
% nz=37;
nt=365;

namef=['N1p';'N3n'];

% fields for BFMcoupler_N1pSurfForcFile and BFMcoupler_N3nSurfForcFile
% N1p=1.1055 x 10^(-8) mmol/m^2/s
% N3n=6.5938 x 10^(-7) mmol/m^2/s
dep=[1.1055e-8 6.5938e-7]; % [mmol/m^2/s]

for var=1:2;
    eval(['fidout_',namef(var,:),...
        '=fopen(''',namef(var,:),'_surface_fluxes_GoT_iNEST_V2.dat'',''w'',''l'');']);
    for t=1:nt;
        sflux=ones(nx,ny)*dep(var);
        eval(['fwrite(fidout_',namef(var,:),',sflux,''float32'');']);
    end
    figure;pcolor(sflux');shading flat;colorbar
end
