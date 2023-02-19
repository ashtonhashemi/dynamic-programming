
classdef BatteryModel
	properties
		fPath
		Prm
		Model
		Capacity
	end
	methods
		% Class Initialize %
		function My= BatteryModel()
			
			% Set Model Type %
			My.Model= 'RC Simple';
			
			% Set Characteristics File %
			My.fPath= 'BatteryLibrary.XLSX';
			
			% Check The Need To Add Data File Path %
			if exist(My.fPath)== 0
				My.fPath= ['../../Parameter/' My.fPath];
			end
			
			% Read Characteristics From Excel File %
			My= My.ReadCharacteristics();
			
			% PlotCharacteristics(My, 'Print');
			% PlotBehavior(My);
			
		end
		function My= ReadCharacteristics(My)
			% Read Data File %
			disp '<< Read Parameters: Battery >>'
			
			% Read Global Data %
			My.Prm.Ns		= 4;
			My.Prm.Np		= 1;
			
			% Read Special Model Data %
			switch My.Model
				case 'RC Afshin'
					
					% Load Temporary Battery Data %
					Bat = load('PMP_Sim_Bat_Li_Ion_7_290');
					
					% Fit Curve To Experimental Data %
					Bat				= BatteryManipulator(Bat);
					My.Prm.Q		= Bat.Q_cap * 3600 *0.05;
					My.Prm.SIn		= Bat.SOC_Ind;
					My.Prm.TIn		= Bat.Tmp_Ind;
					My.Prm.Rc		= Bat.Res_cha;
					My.Prm.Rd		= Bat.Res_dis;
					My.Prm.Voc		= Bat.Vol_cha;
					
				case 'Integrator'
					My.Prm.Eoc	= (My.Prm.Ns)* 12.3688;						% Open Circuit Voltage
					
					My.Prm.Q	= (My.Prm.Np)* 26* 3600* My.Prm.Eoc;		% Charge Capacity [W.Sec]
					My.Capacity	= My.Prm.Q;									% Energy Capacity [Joule]
				case 'RC Simple'
					My.Prm.Ri	= (My.Prm.Ns/My.Prm.Np)* 0.02101;			% Internal Resistance
					My.Prm.Eoc	= (My.Prm.Ns)* 12.3688;						% Open Circuit Voltage
					
					My.Prm.Q	= (My.Prm.Np)* 26* 3600;					% Charge Capacity [A.Sec]
					My.Capacity	= My.Prm.Q* My.Prm.Eoc;						% Energy Capacity [Joule]
				case 'RC Lookup'
					Record			= xlsread(My.fPath, 'VRLA');
					
					My.Prm.SOC	= Record(:,1);								% [0-1]
					My.Prm.Eoc	= (My.Prm.Ns)* Record(:,2);					% [V]
					My.Prm.Rc	= (My.Prm.Ns/My.Prm.Np)* Record(:,3);		% [Ohm]
					My.Prm.Rd	= (My.Prm.Ns/My.Prm.Np)* Record(:,4);		% [Ohm]
					
					My.Prm.Q	= (My.Prm.Np)* 26* 3600;					% Charge Capacity [A.Sec]
					My.Capacity	= My.Prm.Q* My.Prm.Eoc;						% Energy Capacity [Joule]
			end
			
		end
		% Global Function %
		function PlotCharacteristics(My, Flag)
			
			% This Plot Executed Only For Lookup Model %
			if ~strcmp(My.Model, 'RC Lookup'), return; end
			
			figure;
			hAxis(1)= subplot(121);
			hAxis(2)= subplot(122);
			
			PlotVoltage(My, hAxis(1));
			PlotResistance(My, hAxis(2));
			
			if strcmp(Flag, 'Print')
				% Print Characteristics Plot %
				PrintToImage(gcf, [], []);
			end
			
		end
		function PlotOperation(My, T, P, hAxis, Vehicle, Dt)
			
			% This Plot Executed Only For Lookup Model %
			if ~strcmp(My.Model, 'RC Lookup'), return; end
			
			PlotResistance(My, hAxis(1));
			PlotVoltage(My, hAxis(2));
			
			%     % Plot Fuel Cell Operation Point %
			%     [MasOp, EffOp, PowOp]= FuelCellModel(P, Vehicle.FC);
			% 	plot(hAxis(1),     P/ 1E3, EffOp, 'bO', 'MarkerSize', 6);
			% 	plot(hAxis(2), PowOp/ 1E3, MasOp, 'bO', 'MarkerSize', 6);
			%
			%     % plot Equivalent Gasoline / 100 km	%
			% 	PlotEquivalentConsumption(hAxis(3), Vehicle, T, Dt, MasOp);
			
		end
		% Special Function %
		function PlotEfficiency(My, hAxis)
			
			% Prepare Axes For Plot %
			hAxis= PrepareAxis(hAxis, 1);
			
			% Read Battery Data %
			Ri= My.Prm.Ri;
			Eoc= My.Prm.Eoc;
			
			Pt= -5000:1:10000;
			
			% Obtain Current From KVL %
			I= (-Eoc+ (Eoc.^2+ 4* Ri.* Pt).^0.5)./(2* Ri);
			
			iCharge= find(I> 0);
			iDisCharge= find(I<= 0);
			
			Eta(iCharge)= Eoc* I(iCharge)./ Pt(iCharge);
			Eta(iDisCharge)= 1./ (Eoc* I(iDisCharge)./ Pt(iDisCharge));
			
			% Plot Battery Efficiency %
			hold(hAxis, 'on');
			
			plot(hAxis, [Pt(iDisCharge), Pt(iCharge)]/ 1E3, [Eta(iDisCharge), Eta(iCharge)]* 100);
			
			xlabel(hAxis, 'Terminal Power [kW]');
			ylabel(hAxis, 'Efficiency [%]');
			
		end
		function PlotVoltage(My, hAxis)
			
			% Read Battery Data %
			Soc= My.Prm.SOC;
			Eoc= My.Prm.Eoc;
			
			SocTmp= Soc(1):0.001:Soc(end);
			EocTmp= interp1(Soc, Eoc, SocTmp, 'pchip');
			
			% Plot Battery Characteristics %
			hold(hAxis, 'on');
			
			plot(hAxis, SocTmp, EocTmp, 'b-', 'LineWidth', 1);
			plot(hAxis, Soc, Eoc, 'bO', 'MarkerSize', 4);
			
			xlabel(hAxis, 'State Of Charge');
			ylabel(hAxis, 'Open Circuit Voltage [V]');
			set(hAxis, {'Box', 'XGrid', 'YGrid'}, {'On', 'On', 'On'});
			
			set(hAxis, 'YAxisLocation', 'Right');
			
		end
		function PlotResistance(My, hAxis)
			
			% Read Battery Data %
			Soc= My.Prm.SOC;
			Rd = My.Prm.Rd* 1000;
			Rc = My.Prm.Rc* 1000;
			
			SocTmp= Soc(1):0.001:Soc(end);
			RdTmp= interp1(Soc, Rd, SocTmp, 'pchip');
			RcTmp= interp1(Soc, Rc, SocTmp, 'pchip');
			
			% Plot Battery Characteristics %
			hold(hAxis, 'on');
			
			plot(hAxis, SocTmp, RdTmp, 'r-', 'LineWidth', 1);
			hTemp= plot(hAxis, Soc, Rd, 'rO', 'MarkerSize', 4);
			set(get(get(hTemp, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
			
			plot(hAxis, SocTmp, RcTmp, 'b-', 'LineWidth', 1);
			hTemp= plot(hAxis, Soc, Rc, 'bO', 'MarkerSize', 4);
			set(get(get(hTemp, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
			
			xlabel(hAxis, 'State Of Charge');
			ylabel(hAxis, 'Internal Resistance [m\Omega]');
			set(hAxis, {'Box', 'XGrid', 'YGrid'}, {'On', 'On', 'On'});
			
			% Assign Legend %
			hLeg= legend(hAxis, 'Discharge', 'Charge');
			set(hLeg, 'Color', 'None', 'Box', 'Off', 'FontSize', get(0, 'DefaultAxesFontSize')- 1);
			
		end
		function PlotVoltageResistance(My, hAxis)
			
			% 			V= xlsread('Voltage Source.XLS');
			% 			R= xlsread('Internal Resistance.XLS');
			%
			% 			S= 0:10:100;
			% 			V= V(:,2);
			% 			R= R(:,2);
			%
			% 			S1= 0:1:100;
			% 			V1= interp1(S, V, S1, 'spline');
			% 			R1= interp1(S, R, S1, 'spline');
			%
			% 			hold all
			%
			% 			[AX,H1,H2]= plotyy(S, V, S, R);
			% 			[AX,H3,H4]= plotyy(S1, V1, S1, R1);
			%
			% 			set(AX(1), 'XLim', [0,100]);
			% 			set(AX(2), 'XLim', [0,100]);
			%
			% 			set(get(AX(1),'Ylabel'), 'String', 'Open Circuit Voltage [V]', 'FontName', 'Tahoma', 'FontSize', 12);
			% 			set(get(AX(2),'Ylabel'), 'String', 'Internal Resistance [m\Omega]', 'FontName', 'Tahoma', 'FontSize', 12);
			%
			% 			xlabel('State Of Charge [%]', 'FontName', 'Tahoma', 'FontSize', 12);
			%
			% 			set(H1, 'LineStyle','none', 'LineWidth', 2, 'Marker', 'O', 'markersize',8);
			% 			set(H2, 'LineStyle','none', 'LineWidth', 2, 'Marker', 'S', 'markersize',8);
			%
			%
			% 			set(AX(1),'YColor','b')
			% 			set(AX(2),'YColor','r')
			%
			% 			set([H1, H3], 'Color','b');
			% 			set([H2, H4], 'Color','r');
			%
			% 			set(H3, 'LineStyle','-', 'LineWidth', 1);
			% 			set(H4, 'LineStyle','-', 'LineWidth', 1);
			%
			% 			box on;
			% 			grid on;
		end
		function PlotBehavior(My)
			
			% This Plot Executed Only For Lookup Model %
			if ~strcmp(My.Model, 'RC Simple'), return; end
			
			figure;
			hAxis= axes;
			
			
			Imax= +8;
			Imin= -4;
			
			Pmax= +50;
			Pmin= -100;
			
			I=Imin:0.1:Imax;
			E= 12;
			Ri= 2;
			
			Vb= E- Ri* I;
			
			Pb= E* I- Ri* I.^2;
			
			[AX,H1,H2] = plotyy(I, Pb, I, Vb, 'plot');
			linkaxes(AX, 'y');
			
			set(get(AX(1),'Ylabel'), 'String', 'Pb')
			set(get(AX(2),'Ylabel'), 'String', 'Vb')
			
			set(H1,'LineStyle','--')
			set(H2,'LineStyle',':')
			
			axis(AX(1))
			line([Imin, 0], [0, 0], 'LineWidth',3,'LineStyle','--')
			line([0, 0], [Pmin, 0], 'LineWidth',3,'LineStyle','--')
			text(Imin+ 0.5, -4,' Charge Domain','FontSize',10, 'VerticalAlignment', 'Middle')
			line([0, E/(Ri)], [0, 0], 'LineWidth', 3, 'LineStyle', '--', 'Color', 'red')
			line([0, 0], [0, E^2/(4*Ri)], 'LineWidth',3,'LineStyle','--', 'Color', 'red')
			text(1, 4,' DisCharge Domain', 'FontSize',10, 'VerticalAlignment', 'Middle')
			
			hold on
			plot([Imin Imax],[0 0])
			
			plot([0 0],[Pmin Pmax])
			plot(E/(Ri), 0, 'O')
			plot(E/(2*Ri), E^2/(4*Ri), 'O')
			
			axis(AX(2))
			plot(0, E, 'O')
			
			box on
		end
		function CalcPowerRange(My)
			% Calculate Power Range Based On SOC Range %
			
		end
		% Battery Model %
		function X2= SolveSOC(My, X1, P1, Dt)
			
			switch My.Model
				case 'RC Afshin'
					X2= My.ModelRCAfshin(X1, P1, Dt);
				case 'Integrator'
					X2= My.ModelIntegrator(X1, P1, Dt);
				case 'RC Simple'
					X2= My.ModelRC(X1, P1, Dt);
				case 'RC Lookup'
					X2= My.ModelRCLookup(X1, P1, Dt);
			end
			
		end
		function X2= SolveThermal(My, X1, P1, Dt)
			
			% 			% Simple RC Circuit Model %
			% 			Eb= My.Prm.Eoc;
			% 			Qb= My.Prm.Q;
			% 			Ri= My.Prm.Ri;
			%
			% 			% Obtain Current From KVL %
			% 			I= (-Eb+ (Eb^2+ 4* Ri* P1).^0.5)/(2* Ri);
			%
			%
			% 			% Battery Data %
			% 			m  = 10;
			% 			Cp = 795/8;
			% 			hc = 1;
			% 			A  = 1;
			% 			Tamb = 25 + 272.15;
			%
			% 			% State: Temperature %
			% 			Qc = -P1;
			%
			% 			X2= X1+ Dt*( (Ri*I.^2 - hc*A.*(X1 - Tamb) + Qc) /(m * Cp)) ;
			
		end
		function [X2]= ModelIntegrator(My, X1, P1, Dt)
			% Simple Integrator Model %
			Q= My.Prm.Q;
			
			% State: State Of Charge %
			X2= X1+ (Dt/Q)* (P1);
		end
		function [X2, I]= ModelRC(My, X1, P1, Dt)
			% Simple RC Circuit Model %
			Eb= My.Prm.Eoc;
			Qb= My.Prm.Q;
			Ri= My.Prm.Ri;
			
			% Obtain Current From KVL %
			I= (-Eb+ (Eb^2+ 4* Ri* P1).^0.5)/(2* Ri);
			
			% State: State Of Charge %
			X2= X1+ (Dt/Qb)* (I);
		end
		function [X2, I]= ModelRCLookup(My, X1, P1, Dt)
			
			% Lookup Based RC Circuit Model %
			Eb= interp1(My.Prm.SOC, My.Prm.Eoc, X1, 'pchip');
			
			Rc= interp1(My.Prm.SOC, My.Prm.Rc, X1, 'pchip');
			Rd= interp1(My.Prm.SOC, My.Prm.Rd, X1, 'pchip');
			Ri= (P1< 0).* Rc+ (P1>= 0).* Rd;
			
			Qb= My.Prm.Q;
			
			% Obtain Current From KVL %
			I= (-Eb+ (Eb.^2+ 4* Ri.* P1).^0.5)./(2* Ri);
			
			% State: State Of Charge %
			X2= X1+ (Dt/Qb)* (I);
		end
		function [X2, I]= ModelRCAfshin(My, X1, P1, Dt)
			
			% Data Manipulation %
			% 			Pe = U{1};
			% 			Qc = U{2};
			
			
			
			S  = X1{1};
			T  = X1{2};
			
			% 			Pb = prb.Pd(ti) - Pe;
			Pb= P1{1};
			Qc= P1{2};
			[V,R] = My.BattInterp(S, T, Pb);
			
			
			% Battery Data %
			m  = 10;
			Cp = 795/8;
			hc = 1;
			A  = 1;
			Tamb = 25 + 272.15;
			
			% Battery Temp Dynamic %
			I = -(V - sqrt(V.^2 + 4*R.*Pb))./(2*R);
% 			I = (V - sqrt(V.^2 - 4*R.*Pb))./(2*R);
			X2{2} = X1{2} + Dt*( (R.*I.^2 - hc*A.*(T - Tamb) + Qc) /(m * Cp)) ;
			
			% Battery SOC Dynamic %
			X2{1} = X1{1} + Dt*(I./My.Prm.Q);
			
			
		end
		function [V, R] = BattInterp(My, S, T, Pb)
			T  = T - 273.15;
			V  = interp1(My.Prm.SIn, My.Prm.Voc,  S);
			R  = (Pb <= 0).*(interp2(My.Prm.SIn, My.Prm.TIn, My.Prm.Rc',   S, T, 'linear',max(max(My.Prm.Rc))))...
				+ (Pb >  0).*(interp2(My.Prm.SIn, My.Prm.TIn, My.Prm.Rd',   S, T, 'linear',max(max(My.Prm.Rd))));
			
			
		end % Set Battery Model Values at Defined SOC and Temperature %
		
		function X= SolveModelDynamics(U, Pd, Nt)
			% Run Model Based On Input (Control & Disturbance) %
			X0= 0.5;
			X= zeros(1,Nt);
			X(1)= X0;
			T= 1;
			DT= 1;
			Tk= -DT;
			
			for i=1:Nt- 1
				Tk= Tk+ DT;
				X(i+1)= BatteryDynamics(X(i), U(i), Tk, Pd, DT);
			end
		end
	end
end

