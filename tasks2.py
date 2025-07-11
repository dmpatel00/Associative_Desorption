from pathlib import Path
from gpaw import GPAW
import taskblaster as tb
import numpy as np

distances={
        'H':1.5,
        'C':2.2,
        'O':2.2,
        'Cl':2.8,
        'N':2.2,
        'F':1.8,
        }
BEEF_latticeconstants={
        'Ag': 2.124291441342774,
        'Al': 2.014371968654416,
        'Au': 2.125005742282843,
        'Cu': 1.8557793462386443,
        'Fe': 1.4342135039949488,
        'Ir': 1.9417086002991968,
        'Ni': 1.7814747787301133,
        'Pd': 2.01198110070034,
        'Pt': 2.0164604326360447,
        'Rh': 1.9524473768987403,
        'W': 1.5988082893039466,
        'Co':[ 2.515516227989769, 4.079530231796037],
        'Os':[ 2.7589329992487936, 4.357713272155947],
        'Re':[ 2.787459301514976, 4.4905075709326345],
        'Ru':[ 2.731233112069678, 4.311485469812148],
        'Zn':[ 2.7086575158025092, 4.920304199496273],
        }

def buildsurf(dimensions,symbol):
    from ase.build import fcc111,bcc110,hcp0001
    from ase.data import reference_states,atomic_numbers
    crystalstructure=reference_states[atomic_numbers[symbol]]['symmetry']
    lat_const=BEEF_latticeconstants[symbol]
    if crystalstructure=='fcc':
        a=2*lat_const
        surface=fcc111(symbol,dimensions,a=a,orthogonal=True)
    elif crystalstructure=='hcp':
        a=lat_const[0]
        c=lat_const[1]
        surface=hcp0001(symbol,dimensions,a=a,c=c,orthogonal=True)
    elif crystalstructure=='bcc':
        a=2*lat_const
        surface=bcc110(symbol,dimensions,a=latconst,orthogonal=True)
    c = FixAtoms(indices = list(range(24)))
    surface.cell[2][2] = 30
    surface.set_constraint(c)
    surface.center(axis=2)
    return surface

@tb.workflow
class SurfWorkflow:
    surface = tb.var()
    calculator = tb.var()

    @tb.task
    def relax(self):
        return tb.node(
            'minimize_surf',
            surface=self.surface,
            calculator=self.calculator)
def minimize_surf(surface,calculator):
    from gpaw import GPAW, Mixer, MixerFull
    from ase.optimize import BFGS
    from ase.db import connect
    symbol=surface[0].symbol
    adsorbate='clean'
    db=connect(f'{symbol}_{adsorbate}.db')
    if symbol in ['Fe','Co','Ni']:
        spinpol=True
        mixer = MixerFull(beta=0.02,
              nmaxold=5,
              weight=50.0)
    else:
        spinpol=False
        mixer = Mixer(beta=0.02,
              nmaxold=5,
              weight=50.0)
    site='None'
    surface.calc=GPAW(**calculator,spinpol=spinpol,mixer=mixer)
    qn = BFGS(surface,logfile='%s_%s_%s.log'%(symbol,adsorbate,site))
    qn.run(fmax = 0.05)
    db.write(surface,element=symbol,adsorbate=adsorbate,site=site)
    return surface

@tb.workflow
class AdsorbateWorkflow:
    surface = tb.var()
    adsorbate = tb.var()
    calculator = tb.var()

    @tb.task
    def relax(self):
        return tb.node(
            'calc_adsorbate',
            surface=self.surface,
            adsorbate=self.adsorbate,
            calculator=self.calculator)

def calc_adsorbate(surface,adsorbate,calculator):
    from ase.build import add_adsorbate
    from ase.io import write
    from gpaw import GPAW, Mixer, MixerFull
    from ase.optimize import BFGS
    from ase.db import connect
    symbol=surface[0].symbol
    db=connect(f'{symbol}_{adsorbate}.db')
    if symbol in ['Fe','Co','Ni']:
        spinpol=True
        mixer = MixerFull(beta=0.02,
              nmaxold=5,
              weight=50.0)
    else:
        spinpol=False
        mixer = Mixer(beta=0.02,
              nmaxold=5,
              weight=50.0)
    sites=surface.info['adsorbate_info']['sites'].keys()
    min_atoms_vec=[]
    for site in sites:
        surf_ads=surface.copy()
        add_adsorbate(surf_ads,adsorbate,distances[adsobate],site)
        surf_ads.calc=GPAW(**calculator,spinpol=spinpol,mixer=mixer)
        qn = BFGS(surf_ads,logfile='%s_%s_%s.log'%(symbol,adsorbate,site))
        qn.run(fmax = 0.05)
        db.write(surf_ads,element=symbol,adsorbate=adsorbate,site=site)
        min_atoms_vec.append(surf_ads)
    return min_atoms_vec


@tb.workflow
class ParametrizableAdsorbateWorkflow:
    symbol = tb.var()
    calculator = tb.var()

    @tb.task
    def surface(self):
        return tb.node('buildsurface',dimensions=(3,4,3),symbol=self.symbol)

    @tb.subworkflow
    def min_surf(self):
        return SurfWorkflow(surface=self.surface,calculator=self.calculator)
    @tb.subworkflow
    def add_H(self):
        return AdsorbateWorkflow(surface=self.surface,adsorbate='H',calculator=self.calculator)
    @tb.subworkflow
    def add_N(self):
        return AdsorbateWorkflow(surface=self.surface,adsorbate='N',calculator=self.calculator)
    @tb.subworkflow
    def add_O(self):
        return AdsorbateWorkflow(surface=self.surface,adsorbate='O',calculator=self.calculator)
    @tb.subworkflow
    def add_C(self):
        return AdsorbateWorkflow(surface=self.surface,adsorbate='C',calculator=self.calculator)
    @tb.subworkflow
    def add_Cl(self):
        return AdsorbateWorkflow(surface=self.surface,adsorbate='Cl',calculator=self.calculator)
    @tb.subworkflow
    def add_F(self):
        return AdsorbateWorkflow(surface=self.surface,adsorbate='F',calculator=self.calculator)


@tb.dynamical_workflow_generator_task
def parametrize_adsorbate_workflow(calculator):
    material_symbols = [
        'Al', 'Fe','Co','Ni','Cu','Zn','Ru','Rh','Pd','Ag','W','Re','Os','Ir','Pt','Au'
    ]

    for symbol in material_symbols:
        yield f'BEEF-{symbol}', ParametrizableAdsorbateWorkflow(
            symbol=symbol, calculator=calculator)
