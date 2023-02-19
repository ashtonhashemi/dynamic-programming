

% model function %

function [Xj, C, Res] = model(Xi, U, ti, prb, mdl, dem)

% Data Manipulation %
Sr  = U{1};
Sb  = Xi{1};
Sc  = Xi{2};
Pd = dem.P(ti);
[b, c] = modelvar(Sb, Sc, Pd, mdl);
% Value Preparation %
Pb = -Sr.*dem.P(ti);
Pc = -(1-Sr)*dem.P(ti);

if ti==579
disp('here')    ;
end

D1 = b.Voc^2 + 4*b.Rs*Pb<0;
D2 = Sc.^2.*c.Vcp^2 + 4*c.Rs*Pc<0;
    

% Battery SOC Dynamic %
Sbd = (-b.Voc + sqrt(b.Voc^2 + 4*b.Rs*Pb))/(2*b.Rs);
% R = (-Voc + sqrt(Voc^2 + 4*Rs*Pb))/(2*Rs*Q);
Xj{1} = Xi{1} + prb.Ts*(Sbd./b.Q);
Xj{1}(imag(Xj{1})~=0) = 2;
% Capacitor SOC Dynamic %
Scd = -Sc.*(1/(2*c.Rs) + 1/c.Rp) + (sqrt(Sc.^2.*c.Vcp^2 + 4*c.Rs*Pc))./(2*c.Vcp*c.Rs);
Xj{2} = Xi{2} + prb.Ts*(Scd/c.C);
Xj{2}(imag(Xj{2})~=0) = 2;
% Xj{2}(Sc.^2.*c.Vcp^2 + 4*c.Rs*Pc<0)=10;
% Consumption Cost %
Y = 1 - 1./(1 + ((Xj{2} - 0.7)./0.28).^80);
Sbd(imag(Sbd)~=0) = prb.BoundCost;
C = (Sbd.^2) + (prb.BoundPen)*Y ;



% if ti==141
%     dem.P(ti)
%    disp('here') ;
% end
sze = size(U{1});
if sze(1)==1 && sze(2)==1
    Res.X1  = Xi{1};
    Res.X2  = Xi{2};
    Res.Pb  = Pb;
    Res.Pc  = Pc;
    Res.U   = Sr;
    Res.Ib  = Sbd;
    Res.Ic  = Scd;
    % Res.Tm  = Tm;
    % Res.Tg  = Tg;
    % Res.SOC = Xi{1};
else
    Res = NaN;
end


end % Main Model Dynamic %

