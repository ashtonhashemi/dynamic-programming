
function Res = pmp(obj)
UU = 0:0.0005:1;
N  = obj.dem.t(end)/obj.Ts + 1;
X0 = [0.8 0.7];
W  = obj.dem.P;
% P01V = -1000:50:3000;
% P02V = -3000:50:3000;
% obj.P0 = [0 0];
% for II = 1:length(P01V)
%     % 				(II/length(P01V))*100
%     for JJ = 1:length(P02V)
%         obj.P0 = [P01V(II); P02V(JJ)];
% Call pontryagins main calculation core %
[X, P, U] = PontryaginCore(obj, X0, N, UU, W);
% Final(II,JJ) = X(2,end);
% II*JJ
%     end
% %
% end
% figure;
% surf(P01V, P02V, Final');
% xlabel('P01');
% ylabel('P02');



Res.X1 = X(1,:);
Res.X2 = X(2,:);

Res.U = U;
Res.P = P;
Res.Pb = -Res.U.*W;
Res.Pc = -(1-Res.U).*W;
Q   = 14*26*3600;
C   = 2*500;
Res.Ib = fd1(X, U, Res.Pb, Res.Pc)*Q;
Res.Ic = fd2(X, U, Res.Pb, Res.Pc)*C;


end % initilize PMP and run %
function [X, P, U]  = PontryaginCore(obj, X0, N, UU, W)

P(:,1)    = obj.P0;
X(:,1)    = X0;
% U(1,N) = NaN;

if obj.Mode==1
    UU = 1;
end

for i = 1:N
    Pb = -UU*W(i);
    Pc = -(1-UU)*W(i);
    
    UU(imag(fd1(X(:,i), UU, Pb, Pc))~=0) = NaN;
    UU(imag(fd1(X(:,i), UU, Pb, Pc))~=0) = NaN;
    
    [~,I]      = min(Hamilton(obj, UU, X(:,i), P(:,i), W, Pb, Pc));
    U(:,i)     = UU(I);
    Pb = -U(:,i)*W(i);
    Pc = -(1-U(:,i))*W(i);
    X(:,i + 1) = Model (obj, X(:,i), U(:,i), W, Pb, Pc);
    P(:,i + 1) = GradX (obj, X(:,i), U(:,i), P(:,i), W, Pb, Pc);
    
end % Forward movements in time
% close(h);

end % Main PMP routin minimizing hamiltonian
function H    = Hamilton(obj, U, X, P, W, Pb, Pc)

% H = 4*X(1).^2 + U.^2 + P(1).*X(2) + P(2)*U;
Q   = 14*26*3600;
% H = U.^2 + P(1)*fd1(X, U, Pb, Pc) + P(2)*fd2(X, U, Pb, Pc);
y   = obj.L*(1 - 1./(1 + ((X(2) - 0.7)./0.28).^(2*40)));
H = (fd1(X, U, Pb, Pc)*Q).^2 + y + P(1)*fd1(X, U, Pb, Pc) + P(2)*fd2(X, U, Pb, Pc);

end % Return hamiltionian value
function Xj   = Model(obj, X, U, W, Pb, Pc)

Xj =  X + obj.Ts*[fd1(X, U, Pb, Pc); fd2(X, U, Pb, Pc)];

end % State variables dynamic
function Pj   = GradX(obj, X, U, P, W, Pb, Pc)

z = ((X(2) - 0.7)./0.28);
Y = obj.L*(80/0.28)*((        z.^(79)      )./(1 +     (z).^(2*40)      ).^2);
Pj = P - obj.Ts*[0; (P(2)*fdx2(X, U, P, Pb, Pc) + Y)];

end % Co-state variables dynamic

function R = fd1(X, U, Pb, Pc)

Voc = 24;
Rs  = (1E-3)*26.9*2/14;
Q   = 14*26*3600;
R = (-Voc + sqrt(Voc^2 + 4*Rs*Pb))/(2*Rs*Q);


end
function R = fd2(X, U, Pb, Pc)
Rs = (1/2)*2.1*(1E-3);
Rp = (1/2)*3.08*(1E3);
Vcp = 16;
C   = 2*500;
R =(-X(2)*(0.5/Rs + 1/Rp) + sqrt(Vcp^2*X(2)^2 + 4*Rs*Pc)/(2*Rs*Vcp))/(C);



end
function R = fdx1(X, U, P, Pb, Pc)
R = 0;


end
function R = fdx2(X, U, P, Pb, Pc)
Rs = (1/2)*2.1*(1E-3);
Rp = (1/2)*3.08*(1E+3);
Vcp = 16;
C   = 2*500;

R =(-(0.5/Rs + 1/Rp) + Vcp^2*X(2)/(Vcp*2*Rs*sqrt(Vcp^2*X(2)^2 + 4*Rs*Pc)))/C ;


end