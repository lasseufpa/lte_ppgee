

t = poly2trellis(2,[3 1],3);
%u = randi([0 1], 1, 4);
u = [0 1 0 1];
c = convenc(u,t);

mod = comm.RectangularQAMModulator('ModulationOrder', 4, 'BitInput', true);

sym = step(mod, c');

var = 0.7;
stddev = sqrt(var)
noise = randn(length(sym),1);

rSym = sym + stddev * noise;

demod = comm.RectangularQAMDemodulator('ModulationOrder', 4, 'BitOutput', true, ...
    'DecisionMethod', 'Approximate log-likelihood ratio', 'Variance', var);

llr = step(demod, rSym);

app = bcjr2(t, 0.5*ones(1,4), exp(llr)./(1+exp(llr)));

app < 0.5
