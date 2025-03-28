%%%% create masks %%%%%
opengl neverselect
close all
clear all

curdir=pwd;

nx=450;
ny=300;
nz=37;

cd /share/data/squerin/esperimenti/GoT_iNEST/bathymetry

fidin=fopen('./bathydata_GoT_iNEST_V2/hFacC.data','r','l');
for iz=1:1:nz
    [C,count]=fread(fidin,[nx,ny],'float32');
    hFacC_3D(:,:,iz)=C(1:nx,1:ny); %#ok<SAGROW>
    % figure;pcolor(C');shading flat;colorbar
end
mask=hFacC_3D;
surf=mask(:,:,1);
fclose(fidin);

cd(curdir)

qbottom=zeros(nx,ny); % penultima cella d'acqua (sopra la cella d'acqua di fondo)
bottom=zeros(nx,ny); % cella d'acqua di fondo
datum=zeros(nx,ny); % hFacC della cella d'acqua di fondo
for iz=1:nz-1
    for i=1:nx
        for j=1:ny
            if mask(i,j,iz)==1.0 && mask(i,j,iz+1)~=0.0
                qbottom(i,j)=iz;
            end
        end
    end
end

for i=1:nx
    for j=1:ny
        if qbottom(i,j)~=0
            datum(i,j)=hFacC_3D(i,j,qbottom(i,j)+1);
            bottom(i,j)=qbottom(i,j)+1;
        end
    end
end

save('mask_bottom_GoT_iNEST_V2','bottom');
save('mask_qbottom_GoT_iNEST_V2','qbottom');
save('mask_surface_GoT_iNEST_V2','surf');
save('mask_bhFacC_GoT_iNEST_V2','datum');
