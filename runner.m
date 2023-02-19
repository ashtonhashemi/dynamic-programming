function runner(S,~,~)
dc = get(S,'Value');

% Grid Propertios %
grd.X{1}.h = 0.8;
grd.X{1}.l = 0.1;
grd.X{1}.r = 0.01;

grd.X{2}.h = 1;
grd.X{2}.l = 0.4;
grd.X{2}.r = 0.01;

grd.U{1}.h = 1;
grd.U{1}.l = 0;
grd.U{1}.r = 0.01;

grd.U{2}.h = 1;
grd.U{2}.l = 0;
grd.U{2}.r = 1;


% Power Demand Initialization %
[dem, mdl] = drivetrain(dc);

figure;
plot(dem.P);

% Problem Properties %
prb.Ts = 1;
prb.Tf = dem.t(end);
prb.X{1}.t0 = 0.8;
prb.X{1}.f  = 0;
prb.X{2}.t0 = 0.7;
% prb.X{2}.f  = 0.70000; % ECE-15
prb.X{2}.f = 0.7;

prb.N = prb.Tf/prb.Ts +1 ;
prb.FinalCost = [0, 1E9];
prb.BoundCost = 1E12;
prb.BoundPen  = 0; 

% Call DP %
[X, U, ~, ResDP] = dp2dimension(@model, grd, prb, mdl, dem);





% PMP Properties %
% Pontryagins properties %
% prb.Ts  = 1;
% prb.Us  = [0.1; 0.1];
% prb.Uh  = [1 0];
% prb.Ul  = [0 0];
% prb.X0      = [0.8 0.5];
% prb.N = prb.Tf/prb.Ts + 1;
% prb.P01V = [0 79.995999937164140];
% prb.P02V = [0 40.000665222666080];
% 
obj.dem = dem; 
% obj.P0  = [-100, -2259398.943]; % ECE/15
% obj.L   = 1E6; % ECE-15

obj.P0  = [-100, -2513000];
obj.L   = 1E0;

obj.Mode = 0;
obj.Ts  = 1;
Res = pmp(obj);



obj.Mode = 1;
ResE = pmp(obj);

Col = [217 83 25]/255;

t = 0:prb.Ts:prb.Tf;
XLIM = [0 prb.Tf];
close all;
% Graphics %
subplot(4,3,1); box on
hold on
plot(dem.V, 'linewidth',0.5, 'Color', Col);
ylabel('V [km/h]');


subplot(4,3,12); box on
hold on
X = 0.3:0.01:1.1;
Y = 1 - 1./(1 + ((X - 0.7)./0.28).^80);
z = ((X - 0.7)./0.28);
Yp =(80/0.28)*((        z.^(79)      )./(1 +     (z).^(2*40)      ).^2);
[AX,H1,H2] = plotyy(X, Y,  X, Yp);
ylabel('Member Function');
AX(1).XLim = ([0.3 1.1]);
AX(2).XLim = ([0.3 1.1]);

subplot(4,3,4)
hold on
plot(Res.X1, 'linewidth',0.5, 'Color', Col);
plot(ResDP.X1,'--', 'linewidth',0.5, 'Color', 'blue');
% ylim([0.3 0.8]);
% xlabel('Time [s]');
ylabel('S_b');
set(gca,'linewidth',0.5);
box on;
xlim(XLIM);
legend('PMP', 'DP');
title(num2str(Res.X1(end),8));

subplot(4,3,7)
hold on;
plot(Res.X2, 'linewidth',0.5, 'Color', Col);
plot(ResDP.X2,'--', 'linewidth',0.5, 'Color', 'blue');
% xlabel('Time [s]');
ylabel('S_c');
set(gca,'linewidth',0.5);
% ylim([0.4 0.7]);
box on;
xlim(XLIM);
title([num2str(Res.X2(end),8)]);


subplot(4,3,2)
hold on;
plot(Res.U, 'linewidth',0.5,'Color', Col);
plot(ResDP.U,'--', 'linewidth',0.5, 'Color', 'blue');
box on;
% xlabel('Time [s]');
ylabel('r');
xlim(XLIM);

subplot(4,3,3)
hold on
plot(Res.Pb/1000, 'linewidth',0.5, 'Color', Col);
plot(ResDP.Pb/1000,'--', 'linewidth',0.5, 'Color', 'blue');
box on;
% xlabel('Time [s]');
ylabel('P_b_t [kw]');
set(gca,'linewidth',0.5);
xlim(XLIM);
title([num2str(sum(Res.Pb/1000)/3600), ' kwh']);

subplot(4,3,6)
hold on;
plot(Res.Pc/1000,  'linewidth',0.5, 'Color', Col);
plot(ResDP.Pc/1000,'--',  'linewidth',0.5, 'Color', 'blue');
box on;
xlabel('Time [s]');
ylabel('P_c_t [kw]');
set(gca,'linewidth',0.5);
% legend('Dem','Bat','Cap');
box on;
xlim(XLIM);
title([num2str(sum(Res.Pc/1000)/3600), ' kwh']);

subplot(4,3,5)
hold on;
plot(Res.P(1,:),  'linewidth',0.5, 'Color', Col);
plot(Res.P(2,:),'--',  'linewidth',0.5, 'Color', 'blue');
box on;
xlabel('Time [s]');
ylabel('Co-State');
set(gca,'linewidth',0.5);
legend('\lambda_1','\lambda_2');
box on;
xlim(XLIM);
a = Res.X1*(24*14*26*3600) + Res.X2.*((1000*16).^2/1000);
b = (24*14*26*3600) + (1000*16.^2);
c = a/b;
d = a/b - (c(1) - 0.8);
subplot(4,3,10)
hold on
plot(ResE.X1, 'linewidth',0.5, 'Color', Col);
plot(d, 'linewidth',0.5, 'Color', 'black');
% ylim([0.3 0.8]);
% xlabel('Time [s]');
ylabel('S_P_a_c_k');
set(gca,'linewidth',0.5);
box on;
xlim(XLIM);
title(['Single: ',num2str(ResE.X1(end),8), '    Hybrid: ',num2str(d(end),8)]);
legend('Single','Hybrid')

subplot(4,3,9)
hold on
plot(ResE.Pb/1000, 'linewidth',0.5, 'Color', Col);
box on;
% xlabel('Time [s]');
ylabel('P_b_t [kw] (Single)');
set(gca,'linewidth',0.5);
xlim(XLIM);
title([num2str(sum(ResE.Pb/1000)/3600),' kwh']);

ResE.Ib = ResE.Ib(ResE.Ib~=0);
Res.Ib  = Res.Ib(Res.Ib~=0);

Edges = -400:50:100;
subplot(4,3,8)
% plot(ResDP.1Ic);
[N,E] = histcounts(ResE.Ib, Edges);
X = (E(1:end-1) + E(2:end))/2;
m = bar(X, N.*X/3600);
xL = get(gca,'XTickLabel');
xP = get(gca,'XTick');
% Lim = get(gca,'YLim') ;
% histogram(ResDP.Ib);
box on;
xlabel('I_b [A]');
ylabel('I_b [Ah](Single)');
xlim([-400 100]);

subplot(4,3,11)
hold on;
[N,E] = histcounts(Res.Ib, Edges);
X = (E(1:end-1) + E(2:end))/2;
f = bar(X, N.*X/3600);
 set(gca,'XTickLabel', xL);
 set(gca,'XTick', xP);
%  set(gca,'YLim', Lim) ;

% histogram(ResDP.Ib);
box on;
xlabel('I_b [A]');
ylabel('I_b [Ah] (Hybrid)');
xlim([-400 100]);


set(gcf,'PaperUnits','inches')
set(gcf,'Units', 'normalized','PaperType','A4',....
    'PaperOrientation','portrait',...
    'PaperPositionMode','manual');
% ,'PaperPosition',[0 0 3 6]
print(gcf, '-dpng','-r300', 'Result');

ResPMP = Res;
ResElectric = ResE;
save(['Data',dc],'ResPMP', 'ResDP', 'ResElectric');


end



