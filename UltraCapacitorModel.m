
classdef UltraCapacitorModel
	properties
		Prm
		Model
		Capacity
	end
	methods
		% Class Initialize %
		function My= UltraCapacitorModel()
			
			% Set Model Type %
			My.Model= 'Integrator';
			
			% Read Characteristics From Excel File %
			My= My.ReadCharacteristics();
			
		end
		function My= ReadCharacteristics(My)
			% Read Data File %
			disp 'Read Parameters: Ultra Capacitor ...'
			
			% Read Global Data %
			My.Prm.Ns		= 3;
			My.Prm.Np		= 2;
			
			% Read Special Model Data %
			switch My.Model
				case {'Integrator', 'Capacitor'}
					My.Prm.C		= (My.Prm.Np/My.Prm.Ns)* 500;		% Capacitance
					My.Prm.Er		= (My.Prm.Ns)* 16;					% Rated Voltage
					My.Prm.Qc		= 0.5* My.Prm.C* My.Prm.Er^2;		% Charge Capacity [W.s]
					My.Capacity		= My.Prm.Qc;						% Energy Capacity [Joule]
				case 'RRC Circuit'
					My.Prm.Rs		= (My.Prm.Ns/My.Prm.Np)* (2.1/1E3);	% Series Resistance
					My.Prm.Rp		= (My.Prm.Ns/My.Prm.Np)* (3076.9);	% Parallel Resistance

					My.Prm.Rs		= (My.Prm.Ns/My.Prm.Np)* (1.1/1E3);	% Series Resistance
					My.Prm.Rp		= (My.Prm.Ns/My.Prm.Np)* (93076.9);	% Parallel Resistance
					
					My.Prm.C		= (My.Prm.Np/My.Prm.Ns)* 500;		% Capacitance
					My.Prm.Er		= (My.Prm.Ns)* 16;					% Rated Voltage
					
					My.Capacity		= 0.5* My.Prm.C* My.Prm.Er^2;		% Energy Capacity [Joule]
			end
			
		end
		% Global Function %
		function Plot(My)
			
		end
		% Capacitor Model %
		function X2= Solve(My, X1, P1, Dt)
			
			switch My.Model
				case 'Integrator'
					X2= My.ModelIntegrator(X1, P1, Dt);
				case 'Capacitor'
					X2= My.ModelCapacitor(X1, P1, Dt);
				case 'RRC Circuit'
					X2= My.ModelRRC(X1, P1, Dt);
			end
		end
		function X2= ModelIntegrator(My, X1, P1, Dt)
			% Simple Integrator Model %
			Qc= My.Prm.Qc;
			
			% State: State Of Charge %
			X2= X1+ (Dt/Qc)* (P1);
		end
		function X2= ModelCapacitor(My, X1, P1, Dt)
			% Simple Capacitor Model %
			Qc= My.Prm.Qc;
			
			% State: State Of Charge %
 			X2= X1+ (Dt/Qc)* (P1).* (0.5./ X1);
		end
		function X2= ModelRRC(My, X1, P1, Dt)
			% Simple RRC Circuit Model %
			Rs= My.Prm.Rs;
			Rp= My.Prm.Rp;
			C = My.Prm.C;
			Er= My.Prm.Er;
			
			% Calculate Intermediate Parameter %
			A= 1/(2* Rs)+ 1/(Rp);
			B= 1/(2* Rs)* (X1.^2+ 4* Rs* P1/ Er^2).^0.5;
			Sr= (1/C)* (-A* X1+ B);
			
			% State: State Of Charge %
			X2= X1+ Dt* Sr;
		end
	end
end

