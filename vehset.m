

function Veh = vehset
% vehicle mass and aerodynamic data %
Veh.rho   = 1.202;         
Veh.cd    = 0.51;
Veh.af    = 2; 
Veh.mio   = 0.009;
Veh.g     = 9.8;
Veh.m     = 1080;

% Vehicle kinematic Data %
Veh.rt = 0.33;
Veh.RG = [38/11, 43/23, 32/25, 41/39, 35/47];
Veh.rd = 61/17;

end