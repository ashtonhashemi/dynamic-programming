% This function solves two dimensional optimal control problems
% Author: Afshin Pedram
% Date: 3/6/2014
% Version: 1


function  [X, U, Ctg, Res] = dp2dimension(mdl, grd, prb, obj, dem)



% Create State and Control Discretization %
XV{1} = grd.X{1}.l:grd.X{1}.r:grd.X{1}.h;
XV{2} = grd.X{2}.l:grd.X{2}.r:grd.X{2}.h;
UV{1} = grd.U{1}.l:grd.U{1}.r:grd.U{1}.h;
% UV{1} = UV{1}
UV{2} = grd.U{2}.l:grd.U{2}.r:grd.U{2}.h;

% DP backward and forward call
[Cnt, Ctg]  = backwardDP(mdl, XV, UV, prb, obj, dem);
[X, U, Res] = forwardDP(mdl, XV, Cnt, prb, obj, dem);

end  % Dynamic Programming Main Body
function  [Cnt, Ctg]  = backwardDP(mdl, XV, UV, prb, obj, dem)

N = prb.N;
% Initialize Cost-to-go and Control Arrays
Ctg    = zeros(length(XV{1}), length(XV{2}), N + 1);
Cnt{1} = zeros(length(XV{1}), length(XV{2}), N + 1);
Cnt{2} = zeros(length(XV{1}), length(XV{2}), N + 1);
Cnt{1}(:,:,N+1) = 0;
Cnt{2}(:,:,N+1) = 0;
[XD1,XD2] = ndgrid(XV{1},XV{2});
[XX{1}, XX{2}, UU{1}, UU{2}] = ndgrid(XV{1}, XV{2}, UV{1}, UV{2});


% Final State deviation Cost %
Ctg(:,:,N+1) = prb.FinalCost(1)*(XD1-prb.X{1}.f).^2 + prb.FinalCost(2)*(XD2 - prb.X{2}.f).^2 ;


 h = waitbar(0,'Backward DP');
for ti = N:-1:1
    [XXj, Stg, ~] = mdl(XX, UU, ti, prb, obj, dem);
    C = Stg + interpn(XV{1}, XV{2}, Ctg(:, :, ti+1), XXj{1}, XXj{2}, 'linear', prb.BoundCost);
    Shape = [length(XV{1}), length(XV{2}),length(UV{1})*length(UV{2})];
    Crsp = reshape(C,Shape);
    [Ctg(:,:,ti), IndMat] = min(Crsp,[],3);
    Ind{2} = ceil(IndMat/length(UV{1}));
    Ind{1} = IndMat - (Ind{2} - 1)*length(UV{1});
    Cnt{1}(:, :, ti) = UV{1}(Ind{1});
    Cnt{2}(:, :, ti) = UV{2}(Ind{2});
    
    waitbar((N-ti+1)/N,h)
%     clc;
%     Text = sprintf('backward dp: %d %s', floor(100*ti/N), '%');
%     disp(Text);
    
end

delete (h);

end % Backward Dynamic Programming %
function  [X, U, Res]      = forwardDP(mdl, XV, Cnt, prb, obj, dem)

N = prb.N;
X{1}(1) = prb.X{1}.t0;
X{2}(1) = prb.X{2}.t0;
U{1} = zeros(1,N + 1); U{1}(end)=0;
U{2} = zeros(1,N + 1); U{2}(end)=0;
Res.X1(1)  = X{1}(1);
Res.X2(1)  = X{2}(1);
% Res.Wg(1)  = 0;
% Res.Te(1)  = 0;
% Res.Tm(1)  = 0;
% Res.Tg(1)  = 0;
% Res.SOC(1) = X{1}(1);
h = waitbar(0,'Forward DP');
for i = 1:N
    
    U{1}(i) = interp2(XV{1}, XV{2}, Cnt{1}(:,:,i)', X{1}(i), X{2}(i),'linear', NaN);
    U{2}(i) = interp2(XV{1}, XV{2}, Cnt{2}(:,:,i)', X{1}(i), X{2}(i),'linear', NaN);
    
    Xi{1} = X{1}(i);
    Xi{2} = X{2}(i);
    Uf{1} = U{1}(i);
    Uf{2} = U{2}(i);
    [Xj, ~, MyRes] = mdl(Xi, Uf, i, prb, obj, dem);
    X{1}(i + 1) =  Xj{1};
    X{2}(i + 1) =  Xj{2};
    
    Res.X1(i+1)  = MyRes.X1;
    Res.X2(i+1)  = MyRes.X2;
    Res.Pb(i)  = MyRes.Pb;
    Res.Pc(i)  = MyRes.Pc;
    Res.U(i)  = MyRes.U;
    Res.Ib(i)  = MyRes.Ib;
    Res.Ic(i)  = MyRes.Ic;
%     Res.Tm(i+1)  = MyRes.Tm;
%     Res.Tg(i+1)  = MyRes.Tg;
%     Res.SOC(i+1) = MyRes.SOC;
    
    waitbar(i/N,h)
%     clc;
%     Text = sprintf('forward dp: %d %s', floor(100*i/N), '%');
%     disp(Text);
    
    
end
delete(h);

X = [X{1}; X{2}];
U = [U{1}; U{2}];

end % Forward Dynamic Programming %

