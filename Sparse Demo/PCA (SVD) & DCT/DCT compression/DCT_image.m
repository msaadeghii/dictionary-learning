clear all
close all
clc

%
% Read images
I1 = imread('lena.png');I1 = im2double(I1);
I2 = imread('barbara.png');I2 = im2double(I2);
I3 = imread('boat.png');I3 = im2double(I3);
I4 = imread('house.png');I4 = im2double(I4);

[m n]=size(I1);
ne=m*n;

% Apply DCT
dcti1=dct2(I1);
dcti2=dct2(I2);
dcti3=dct2(I3);
dcti4=dct2(I4);

% Set the threshold
Tr=0.05;


% Set small coefficients of DCT to zero

% Percentage of the coefficients becoming zero
ind1=abs(dcti1)<Tr;
n1=sum(ind1(:))/ne;

dcti1(ind1) = 0;
I1new = idct2(dcti1);
figure, subplot(121); imshow(I1);title('Original');
subplot(122); imshow(I1new);title('Compressed');

% Percentage of the coefficients becoming zero
ind2=abs(dcti2)<Tr;
n2=sum(ind2(:))/ne;

dcti2(ind2) = 0;
I2new = idct2(dcti2);
figure, subplot(121); imshow(I2);title('Original');
subplot(122); imshow(I2new);title('Compressed');

% Percentage of the coefficients becoming zero
ind3=abs(dcti3)<Tr;
n3=sum(ind3(:))/ne;

dcti3(ind3) = 0;
I3new = idct2(dcti3);
figure, subplot(121); imshow(I3);title('Original');
subplot(122); imshow(I3new);title('Compressed');

% Percentage of the coefficients becoming zero
ind4=abs(dcti4)<Tr;
n4=sum(ind4(:))/ne;

dcti4(ind4) = 0;
I4new = idct2(dcti4);
figure, subplot(121); imshow(I4);title('Original');
subplot(122); imshow(I4new);title('Compressed');