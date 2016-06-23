import ROOT
import array

ROOT.gSystem.Load('libRooFit.so')
ROOT.gSystem.Load('libFitTools.so')

fitter = ROOT.EfficiencyFitter.singleton()

tree = ROOT.TTree('test', 'test')
x = array.array('d', [0.])
p = array.array('I', [0])
tree.Branch('x', x, 'x/D')
tree.Branch('p', p, 'p/i')

data = []
data += [(0., 1)] * 4
data += [(0., 0)] * 1
data += [(0.5, 1)] * 3
data += [(0.5, 0)] * 2
data += [(1., 1)] * 2
data += [(1., 0)] * 3

for datum in data:
    x[0] = datum[0]
    p[0] = datum[1]
    tree.Fill()

xvar = ROOT.RooRealVar('x', 'x', -2., 2.)
intercept = ROOT.RooRealVar('intercept', 'intercept', 0.2, -10., 10.)
slope = ROOT.RooRealVar('slope', 'slope', 0.5, -1., 1.)
line = ROOT.RooPolyVar('line', 'line', xvar, ROOT.RooArgList(intercept, slope))

fitter.setTarget(tree, 'p', 'x')
fitter.setLikelihood(line, ROOT.RooArgList(intercept, slope))

fitter.fit(1)
