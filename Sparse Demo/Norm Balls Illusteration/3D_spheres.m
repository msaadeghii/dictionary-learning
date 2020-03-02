% Figure 1.2 of Elad's Book. This code is from M. Elad's homepage
% =========================================
% This program demonstrates in 3D the intersection of various Lp 
% spheres with a hyperplane

function []=Chapter_01_3D_spheres() 

% L2 ball intersection with a hyperplane

res=0.01;
p=2;
[x,y,z]=meshgrid(-1-res:res:1+res, ...
   -1-res:res:1+res, -1-res:res:1+res);
w=(abs(x).^p+abs(y).^p+abs(z).^p).^(1/p);
sphereP=isosurface(x,y,z,1.19*w,1);

figure(1); clf; 
count=renderpatch(sphereP,1);
h=patch([-1.5 -1.5 1.5 1.5],[1.5 -1.5 -1.5 1],...
              1.2+0.65*[-1.5 -1.5 1.5 1.5 ]-0.2,[0.8 0 1]);
set(h,'FaceAlpha',0.5); 
hold on;
h=line([-15 15],[0 0],[0,0]); set(h,'LineWidth',2,'LineStyle','-'); 
h=line([0 0],[-15 15],[0,0]); set(h,'LineWidth',2,'LineStyle','-'); 
h=line([0 0],[0,0],[-15 15]); set(h,'LineWidth',2,'LineStyle','-'); 
axis off;
grid on
daspect([1 1 1])
light('position',[10,-10,10])
set(gca,'projection','perspective')
set(gca,'CameraViewAngle',4)
set(gcf,'color', [1 1 1])
view(-4,10)
axis([   -1.5000    1.5000   -2.0000    2.0000   -1.0000    1.0000]);
drawnow;
% print -depsc2 Chapter_01_SpherePlan.eps

% L1.5 ball intersection with a hyperplane

res=0.01;
p=1.5;
[x,y,z]=meshgrid(-1-res:res:1+res, ...
   -1-res:res:1+res, -1-res:res:1+res);
w=(abs(x).^p+abs(y).^p+abs(z).^p).^(1/p);
sphereP=isosurface(x,y,z,1.082*w,1);

figure(1); clf; 
count=renderpatch(sphereP,1);
h=patch([-1.5 -1.5 1.5 1.5],[1.5 -1.5 -1.5 1],...
              1.2+0.65*[-1.5 -1.5 1.5 1.5 ]-0.2,[0.8 0 1]);
set(h,'FaceAlpha',0.5); 
hold on;
h=line([-15 15],[0 0],[0,0]); set(h,'LineWidth',2,'LineStyle','-'); 
h=line([0 0],[-15 15],[0,0]); set(h,'LineWidth',2,'LineStyle','-'); 
h=line([0 0],[0,0],[-15 15]); set(h,'LineWidth',2,'LineStyle','-'); 
axis off;
grid on
daspect([1 1 1])
light('position',[10,-10,10])
set(gca,'projection','perspective')
set(gca,'CameraViewAngle',4)
set(gcf,'color', [1 1 1])
view(-4,10)
axis([   -1.5000    1.5000   -2.0000    2.0000   -1.0000    1.0000]);
drawnow;
% print -depsc2 Chapter_01_L15Plan.eps

% L1 ball intersection with a hyperplane

res=0.01;
p=1;
[x,y,z]=meshgrid(-1-res:res:1+res, ...
   -1-res:res:1+res, -1-res:res:1+res);
w=(abs(x).^p+abs(y).^p+abs(z).^p).^(1/p);
sphereP=isosurface(x,y,z,0.99*w,1);

figure(1); clf; 
count=renderpatch(sphereP,1);
h=patch([-1.5 -1.5 1.5 1.5],[1.5 -1.5 -1.5 1],...
              1.2+0.65*[-1.5 -1.5 1.5 1.5 ]-0.2,[0.8 0 1]);
set(h,'FaceAlpha',0.5); 
hold on;
h=line([-15 15],[0 0],[0,0]); set(h,'LineWidth',2,'LineStyle','-'); 
h=line([0 0],[-15 15],[0,0]); set(h,'LineWidth',2,'LineStyle','-'); 
h=line([0 0],[0,0],[-15 15]); set(h,'LineWidth',2,'LineStyle','-'); 
axis off;
grid on
daspect([1 1 1])
light('position',[10,-10,10])
set(gca,'projection','perspective')
set(gca,'CameraViewAngle',4)
set(gcf,'color', [1 1 1])
view(-4,10)
axis([   -1.5000    1.5000   -2.0000    2.0000   -1.0000    1.0000]);
drawnow;
% print -depsc2 Chapter_01_L1Plan.eps

% Lp ball for p=0.5
res=0.01;
p=0.7;
[x,y,z]=meshgrid(-1-res:res:1+res, ...
   -1-res:res:1+res, -1-res:res:1+res);
w=(abs(x).^p+abs(y).^p+abs(z).^p).^(1/p);
sphereP=isosurface(x,y,z,0.93*w,1);

figure(1); clf; 
count=renderpatch(sphereP,1);
h=patch([-1.5 -1.5 1.5 1.5],[1.5 -1.5 -1.5 1],1.2+0.65*[-1.5 -1.5 1.5 1.5 ]-0.2,[0.8 0 1]);
set(h,'FaceAlpha',0.5); 
hold on;
h=line([-15 15],[0 0],[0,0]); set(h,'LineWidth',2,'LineStyle','-'); 
h=line([0 0],[-15 15],[0,0]); set(h,'LineWidth',2,'LineStyle','-'); 
h=line([0 0],[0,0],[-15 15]); set(h,'LineWidth',2,'LineStyle','-'); 
axis off;
grid on
daspect([1 1 1])
light('position',[10,-10,10])
set(gca,'projection','perspective')
set(gca,'CameraViewAngle',4)
set(gcf,'color', [1 1 1])
view(-4,10)
axis([   -1.5000    1.5000   -2.0000    2.0000   -1.0000    1.0000]);
drawnow;
% print -depsc2 Chapter_01_LpPlan.eps

%================================================

function count = renderpatch(objIn,alpha)
%hierarchical render function for structs and cell arrays
%Takes either a cell array or a single struct as input.
%For each struct, can set:
%  facecolor: default=cyan
%  edgecolor: default='none'
%  ambientstrength: default=.6
%  specularstrength: default=.2
%  specularexponent: default=10
%  facelighting: default='phong'
%  diffusestrength: default=.5
%  visible: default='on'

if (iscell(objIn)) %a list of structs
   
   for i=1:length(objIn)

      obj=patch(objIn{i});
      
      if (isfield(objIn{i},'facecolor'))
         fcolor=objIn{i}.facecolor;
      else
         fcolor=[0,1,1];
      end
      
      if (isfield(objIn{i},'edgecolor'))
         ecolor=objIn{i}.edgecolor;
      else
         ecolor='none';
      end
      
      if (isfield(objIn{i},'ambientstrength'))
         ambstr=objIn{i}.ambientstrength;
      else
         ambstr=.6;
      end
      
      if (isfield(objIn{i},'specularstrength'))
         spcstr=objIn{i}.specularstrength;
      else
         spcstr=.2;
      end
      
      if (isfield(objIn{i},'specularexponent'))
         spcexp=objIn{i}.specularexponent;
      else
         spcexp=10;
      end
      
      if (isfield(objIn{i},'facelighting'))
         facelight=objIn{i}.facelighting;
      else
         facelight='phong';
      end
      
      if (isfield(objIn{i},'diffusestrength'))
         difstr=objIn{i}.diffusestrength;
      else
         difstr=.5;
      end
  
      if (isfield(objIn{i},'visible'))
         vis=objIn{i}.visible;
      else
         vis='on';
      end


      set(obj, 'FaceColor', fcolor, ...
              'FaceAlpha',alpha, ...
               'EdgeColor', ecolor, ...
               'AmbientStrength',ambstr,...
               'SpecularStrength',spcstr, ...
               'SpecularExponent', spcexp, ...
               'FaceLighting', facelight, ...
               'DiffuseStrength', difstr, ...
               'Visible',vis);
   end  
   count=i;
   
 elseif (isstruct(objIn)) %must be a single struct   
    obj=patch(objIn);
    
    if (isfield(objIn,'facecolor'))
         fcolor=objIn.facecolor;
      else
         fcolor=[0,1,1];
      end
      
      if (isfield(objIn,'edgecolor'))
         ecolor=objIn.edgecolor;
      else
         ecolor='none';
      end
      
      if (isfield(objIn,'ambientstrength'))
         ambstr=objIn.ambientstrength;
      else
         ambstr=.6;
      end
      
      if (isfield(objIn,'specularstrength'))
         spcstr=objIn.specularstrength;
      else
         spcstr=.2;
      end
      
      if (isfield(objIn,'specularexponent'))
         spcexp=objIn.specularexponent;
      else
         spcexp=10;
      end
      
      if (isfield(objIn,'facelighting'))
         facelight=objIn.facelighting;
      else
         facelight='phong';
      end
      
      if (isfield(objIn,'diffusestrength'))
         difstr=objIn.diffusestrength;
      else
         difstr=.5;
      end
  
      if (isfield(objIn,'visible'))
         vis=objIn.visible;
      else
         vis='on';
      end

      set(obj, 'FaceColor', fcolor, ...
                        'FaceAlpha',alpha, ...
               'EdgeColor', ecolor, ...
               'AmbientStrength',ambstr,...
               'SpecularStrength',spcstr, ...
               'SpecularExponent', spcexp, ...
               'FaceLighting', facelight, ...
               'DiffuseStrength', difstr, ...
               'Visible',vis);
      count=1;   
 end %if   
