clear all
close all

dataFile='/share/data/vdibiagio/MITgcm_BFM/lightExt/xESP_forc_days1to1920.bin';

nx=570;
ny=264;

strdirGrid='/share/scratch/vdibiagio/MITgcmBFM/CLR_365days/results_merged/';
fhfacC = fopen([strdirGrid,'Computational_hFacC'], 'r','b');
MASKhFac=fread(fhfacC,[nx,ny],'float32');
MASKhFac(MASKhFac~=1)=0;
Kext=zeros(1921,nx,ny);
fKext = fopen(dataFile, 'r','b');
for idays=1:1921
    A= fread(fKext,[nx,ny],'float32');
    Kext(idays,:,:)=A.*MASKhFac;
end

MASKNan=MASKhFac;
MASKNan(MASKNan==0)=NaN;
ANNUAL=mean(Kext(93:93+364,:,:),1);
ANNUAL=squeeze(ANNUAL).*MASKNan;

close all
B=(ANNUAL);
figure;pcolor(B');shading flat;colorbar

% ora interpolo sulla griglia GoT
ycOrigin=45.4;
xcOrigin=13.22;
delX=1/768;
delY=1/768;
nxGoT=450;
nyGoT=300;

mlat=ycOrigin:delY:ycOrigin+(nyGoT-1)*delY;
mlon=xcOrigin:delX:xcOrigin+(nxGoT-1)*delX;
[CLAT,CLON]=meshgrid(mlat,mlon);

fhLON = fopen([strdirGrid,'Computational_XC'], 'r','b');
LON=fread(fhLON,[nx,ny],'float32');
fhLAT = fopen([strdirGrid,'Computational_YC'], 'r','b');
LAT=fread(fhLAT,[nx,ny],'float32');

LATnan=LAT.*MASKNan;
LONnan=LON.*MASKNan;

% "salto" molto grande in risoluzione da Med a GoT_iNEST...
varint=griddata(LON(~isnan(MASKNan)),LAT(~isnan(MASKNan)),...
    ANNUAL(~isnan(MASKNan)),CLON,CLAT,'linear'); % nearest è troppo "a blocchi"...

close all
figure;pcolor(CLAT');shading flat;colorbar;
figure;pcolor(CLON');shading flat;colorbar;
figure;pcolor(LATnan');shading flat;colorbar;
figure;pcolor(LONnan');shading flat;colorbar;
figure;pcolor(varint');shading flat;colorbar;
figure;pcolor(ANNUAL');shading flat;colorbar;
caxis([0.06 0.22]);

maskFileGoT='/share/data/squerin/esperimenti/GoT_iNEST/bathymetry/bathydata_GoT_iNEST_V2/hFacC.data';
fr2=fopen(maskFileGoT,'r','l');
maskGoT=fread(fr2,[nxGoT,nyGoT],'float32');
figure;pcolor((varint.*maskGoT)');shading flat;colorbar;
caxis([0.06 0.22])

fw=fopen('Kext_GoT_iNEST_365.dat','w','l');
for timestep=1:365
    fwrite(fw,varint,'float32');
end

fclose('all')
