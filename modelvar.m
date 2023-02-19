
function [b, c] = modelvar(Sb, Sc, Pd, mdl)
c.Rs = (1/2)*2.1*(1E-3);
c.Rp = (1/2)*3.08*(1E+3);
c.Vcp = 16;
c.C   = 2*500;


b.Voc = 24;
b.Rs  = (1e-3)*26.9*2/14;
b.Q   = 14*26*3600;
end