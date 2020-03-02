clear all
close all
clc

addpath(genpath(cd))

im=imread('lena.png');
im=imresize(im,[128 128]);

im=double(im);
[M N]=size(im);
bb=8;
slidingDis=1;

[y,idx] = my_im2col(im,[bb,bb],slidingDis);

ym=mean(y,2);
y=y-repmat(ym,[1,size(y,2)]);

%%

% Dimension Reduction

no_dim=20;
[mappedX, mapping] = compute_mapping(y','PCA',no_dim);

options.ReducedDim=64;
[U S]=myPCA(y',options);
plot(S,'linewidth',2.5);title('Singular-valuse','fontsize',15);
% Reconstruction
yh = reconstruct_data(mappedX, mapping);
yh=yh';

PSNR=20*log10(norm(y,'fro')/norm(y-yh,'fro'));


%% Reconstruct the image

count = 1;
Weight = zeros(M,N);
IMout = zeros(M,N);
blocks=yh+repmat(ym,[1,size(y,2)]);

[rows,cols] = ind2sub(size(im)-bb+1,idx);
for i  = 1:length(cols)
    col = cols(i); row = rows(i);        
    block =reshape(blocks(:,count),[bb,bb]);
    IMout(row:row+bb-1,col:col+bb-1)=IMout(row:row+bb-1,col:col+bb-1)+block;
    Weight(row:row+bb-1,col:col+bb-1)=Weight(row:row+bb-1,col:col+bb-1)+ones(bb);
    count = count+1;
end;


im_comp = IMout./Weight;
ssim=ssim_index(im,im_comp);
psnr=20*log10(norm(im,'fro')/norm(im-im_comp,'fro'));
figure;
subplot(121);imshow(uint8(im));title('Original','fontsize',15);
subplot(122);imshow(uint8(im_comp));title(['PSNR: ',num2str(psnr),', SSIM: ',num2str(ssim)],'fontsize',15);