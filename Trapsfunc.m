x = 0:0.001:1.2;
y = 1 - 1./(1 + ((x - 0.7)./0.28).^(2*40));

plotyy(x,y,x,Y);