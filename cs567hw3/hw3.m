%%
%1

clear;
file_path = '~/Downloads/AMIP_II_CD/LessonData/7_SpatialTransforms';
img=double(imread(sprintf('%s/PorkyPig.jpg', file_path)));

rotimg=zeros(300,300);
angle=-20;
rotmat=zeros(2,2);
rotmat(1,1)=cos(angle*pi/180);
rotmat(1,2)=-sin(angle*pi/180);
rotmat(2,1)=sin(angle*pi/180);
rotmat(2,2)=cos(angle*pi/180);

rotmat_inv=inv(rotmat);
rotimg=zeros(300,300);
oldpix=zeros(2,1);
newpix=zeros(2,1);

shift=zeros(2,1);
shift(1,1)= 150;
shift(2,1)= 150;
rho=0;

for i=1:300
    for j=1:300
        newpix(1,1)=i;
        newpix(2,1)=j;
	    oldpix=round(rotmat_inv*(newpix-shift)+shift);
        if (oldpix(1,1) > 0) & (oldpix(1,1) < 300) & (oldpix(2,1) > 0) & (oldpix(2,1) < 300)
            rho=img(oldpix(1,1),oldpix(2,1));
        end
	rotimg(i,j)=rho;
    end
end

rotbilinear=zeros(300,300);
rho = 0;
for i=1:300
    for j=1:300
        newpix(1,1)=i;
        newpix(2,1)=j;
	    oldpix=rotmat_inv*(newpix-shift)+shift;
        diff = rem(oldpix,1);
        if ((floor(oldpix(1,1)) > 0) & (floor(oldpix(1,1)) < 300) & (floor(oldpix(2,1)) > 0) & (floor(oldpix(2,1)) < 300))
            xvec = [1-diff(1,1);diff(1,1)]';
            ker = [img(floor(oldpix(1,1))+1,floor(oldpix(2,1))),img(floor(oldpix(1,1))+1,floor(oldpix(2,1))+1);img(floor(oldpix(1,1)),floor(oldpix(2,1))),img(floor(oldpix(1,1)),floor(oldpix(2,1))+1)];
            yvec = [1-diff(2,1);diff(2,1)];
            rho=xvec*ker*yvec;
        end
	rotbilinear(i,j)=rho;
    end
end


figure
colormap(gray)
subplot(1, 2, 1)
imagesc(rotimg)
title('Nearest Neighbour')
subplot(1,2,2)
imagesc(rotbilinear)
title('Bilinear Interpolation')


%%
%2
clear;
file_path = '~/Downloads/AMIP_II_CD/LessonData/9_Registration';
imgT1= double(imread(sprintf('%s/T1.jpg', file_path)));
imgCT= double(imread(sprintf('%s/CT.jpg', file_path)));
size(imgT1)
size(imgCT)

% H for T1
jhist=zeros(1,256);
for i=1:388
    for j=1:388
        jhist(1,imgT1(i,j)+1) = jhist(1,imgT1(i,j)+1) + 1;
    end
end

HT1 = 0;
for i=1:256
    if jhist(1,i)~=0
        p = jhist(1,i)./(size(imgT1,1)*size(imgT1,2));
        HT1 = HT1 + -1*p*log(p);
    end
end
HT1

% H for CT
jhist=zeros(1,256);
for i=1:388
    for j=1:388
        jhist(1,imgCT(i,j)+1) = jhist(1,imgCT(i,j)+1) + 1;
    end
end

HCT = 0;
for i=1:256
    if jhist(1,i)~=0
        p = jhist(1,i)./(size(imgCT,1)*size(imgCT,2));
        HCT = HCT + -1*p*log(p);
    end
end
HCT

jhist=zeros(1,256);
figure
for a=-9:10
    rotimg=zeros(388,388);
    angle=a*10;
    rotmat=zeros(2,2);
    rotmat(1,1)=cos(angle*pi/180);
    rotmat(1,2)=-sin(angle*pi/180);
    rotmat(2,1)=sin(angle*pi/180);
    rotmat(2,2)=cos(angle*pi/180);
    oldpix=zeros(2,1);
    newpix=zeros(2,1);
    invrotmat=transpose(rotmat);
    shift=transpose([194,194]);
    
    for i=1:388
        for j=1:388
            newpix(1,1)=i;
            newpix(2,1)=j;
            oldpix=round((invrotmat*(newpix-shift))+shift);
            if (oldpix(1,1) > 0) & (oldpix(1,1) < 388) & (oldpix(2,1) > 0) & (oldpix(2,1) < 388)
                rotimg(i,j)=imgT1(oldpix(1,1),oldpix(2,1));
            end
        end
    end
    
    jhist=zeros(256,256);
    for i=1:388
        for j=1:388
            rho1=imgCT(i,j)+1;
            rho2=rotimg(i,j)+1;
            jhist(rho1,rho2)=jhist(rho1,rho2)+1;
        end
    end
   Hjoint = 0;
    for i=1:256
        for j=1:256
            if jhist(i,j)~=0
                p = jhist(i,j)./(size(imgCT,1)*size(imgCT,2));
                Hjoint = Hjoint + -1*p*log(p);
            end
        end
    end
    mutinfo = (HCT+HT1)./Hjoint;
    m = 4;
    s = 4;
    jhist = 63*1./(1+exp(-1*(jhist-m)/s));
    subplot(4, 5, a+10)
    imagesc(jhist)
    colorbar
    title(sprintf('Joint Entropy = %.4g \n Mutual Info = %.4g', Hjoint,mutinfo));
end

