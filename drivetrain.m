    

function [dem, mdl]= drivetrain(drivIndex)

% read driving cycle %
Veh = vehset;
dri = ReadDrivingCycle(drivIndex);
mdl = ReadNominal;

% calculating demaned power and force %
dem.A     =  diff(dri.v)./diff(dri.t);
Finertial =  Veh.m*(diff(dri.v)./diff(dri.t));
Faero     =  0.5*Veh.rho*Veh.cd*Veh.af*((dri.v(1:end-1) + dri.v(2:end))/2).^2 ;
Ffri      =  Veh.m*Veh.g*Veh.mio*ones(size(Faero));
Fgrade    =  Veh.m*Veh.g*sin(dri.e); 
dem.t     =  dri.t;
dem.V     =  dri.v;
dem.F     =  Finertial + Faero + Ffri + Fgrade;
dem.F     =  (dem.F>0).*dem.F + 0.2*(dem.F<0).*dem.F;

% Tire %
dem.W     = dem.V./Veh.rt;
dem.T     = dem.F.*Veh.rt;

% Differentil %
dem.W     = dem.W.*Veh.rd;
dem.T     = dem.T./Veh.rd;
dem.T = [dem.T 0];
% Gear Box %
dem.V = dem.V*3.6;
rg = Veh.RG(1)*(dem.V<62) + Veh.RG(2)*(dem.V>=62 & dem.V<82) + ...
         Veh.RG(3)*(dem.V>=82 & dem.V<109) + Veh.RG(4)*(dem.V>=109 & dem.V<143) +...
         Veh.RG(5)*(dem.V>=143);
dem.W     = dem.W.*rg;
dem.T     = dem.T./rg;


% Motor %
eff = motoeff(dem.T, dem.W);

% To electric line %
% dem.P     =  eff.*[dem.T.*(dem.W(1:end-1) + dem.W(2:end))/2 0];
dem.P     =  eff.*dem.T.*(dem.W);

% Apply Shifting Strategy and Driveline Efficiency %

end

function dri = ReadDrivingCycle(Index)
drilist = {'US06-HWY','EUDC','NEDC','FTP-75','HWFET',...
           'ECE-15','FTP-HWY','UDDS','10-Mode','15-Mode', 'HWFET-MTN','WLTC1','WLTC2','WLTC3'};  
data    = load(['DrivingCycle\' drilist{Index}]);
dri.t   = 0:1:max(data.t);
dri.v   = interp1(data.t,data.v,dri.t);
if Index == 11 
dri.e   = 10*interp1(data.t,data.e,dri.t(1:end-1));
else
dri.e   = zeros(1, length(dri.t)-1);    
end


end

function obj = ReadNominal

% obj.eng = advisorfile('FC_PRIUS_JPN','Engine');
% obj.mot = advisorfile('MC_PRIUS_JPN','Motor');
% obj.gen = advisorfile('GC_PRIUS_JPN','Gene');
obj = 0;
end

function eff = motoeff(T,W)
data = load('OC_SIM_DB_Mot_pm_50_78_premag');
Mot = data.Mot;

E = interp2(Mot.sp(13:26), Mot.tq(:), 1./Mot.eff(:,13:26), W, abs(T))*0.95*0.98;
eff = E.*(T>=0) + (1./E).*(T<0);


figure;
hold on;
plot(W, T, '*');
plot(Mot.sp_full(15:29), Mot.tq_max(15:29),'--','Color','red');
[a, b] = contour(Mot.sp(14:26), Mot.tq(12:22), 1./Mot.eff(12:22,14:26),[0.94 0.93  0.92 0.9 0.85 0.8 0.75]);
clabel(a,b);

end










