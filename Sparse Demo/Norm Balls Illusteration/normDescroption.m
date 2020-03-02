% Figure 1.3  of Elad's Book. This code is from M. Elad's homepage
% =========================================
% Showing the L0 norm effect to count non-zeros

figure(1); clf; 
x=-2:0.01:2;

p=[2 1 0.5 0.1];
for k=1:1:4,
    h=plot(x,abs(x).^p(k),'k'); 
    set(h,'LineWidth',2);
    hold on;
end;
axis([-1.5 1.5 0 2]); 
set(gca,'FontSize',18);
disp('Position the text p=2 on the figure');
h=gtext('p=2'); 
set(h,'FontSize',18);
disp('Position the text p=1 on the figure');
h=gtext('p=1'); 
set(h,'FontSize',18);
disp('Position the text p=0.5 on the figure');
h=gtext('p=0.5'); 
set(h,'FontSize',18);
disp('Position the text p=0.1 on the figure');
h=gtext('p=0.1'); 
set(h,'FontSize',18);

% print -depsc2 Chapter_01_normDescroption.eps