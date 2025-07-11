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

def getsurface(symbol):
    from ase.db import connect
    db=connect('/home/cat/dmapa/gpaw/dissads/adsorbates.db')
    surface=db.get(element=symbol,adsorbate='clean')
    return surface
ads_map={
        'fcc':{
            'H2':[16,20,36],
            'N2':[16,20,36],
            'O2':[16,20,36],
            'CO':[16,20,36],
            'CN':[16,20,36],
            'F2':[16,20,36],
            'Cl2':[16,20,36],
            }
        }

@tb.workflow
class EndStateWorkflow:
    surface = tb.var()
    molecule = tb.var()
    calculator = tb.var()

    @tb.task
    def initial(self):
        return tb.node(
            'initial_state',
            surface=self.surface,
            molecule=self.molecule,
            calculator=self.calculator)

    @tb.task
    def final(self):
        return tb.node(
            'final_state',
            surface=self.surface,
            molecule=self.molecule,
            calculator=self.calculator)

def initial_state(surface,molecule,calculator):
    from gpaw import GPAW, Mixer, MixerFull
    from ase.optimize import BFGS
    from ase.build import molecule
    from ase.io import write
    from ase.data import reference_states,atomic_numbers
    symbol=surface[0].symbol
    crystalstructure=reference_states[atomic_numbers[symbol]]['symmetry']
    if symbol in ['Fe','Co','Ni']:
        spinpol=True
        mixer = MixerFull(beta=0.02,
              nmaxold=5,
              weight=100.0)
        for atom in surface:
            atom.magmom=1.4
    else:
        spinpol=False
        mixer = Mixer(beta=0.05,
              nmaxold=5,
              weight=100.0)
    adsorbate_indices=ads_map[crystalstructure][molecule]
    surface.append(molecule[0].symbol)
    surface.positions[-1]=surface.positions[adsorbate_indices[0]]+(0,0,distances[molecule[0].symbol])
    surface.append(molecule[1].symbol)
    surface.positions[-1]=surface.positions[adsorbate_indices[1]]+(0,0,distances[molecule[1].symbol])
    surface.calc=GPAW(**calculator,spinpol=spinpol,mixer=mixer)
    qn = BFGS(surface,logfile='%s_%s.log'%(symbol,molecule))
    qn.run(fmax = 0.05)
    write(f'{symbol}_{molecule}_initial.traj'),
    return surface

def final_state(surface,molecule,calculator):
    from gpaw import GPAW, Mixer, MixerFull
    from ase.optimize import BFGS
    from ase.build import molecule as buildmol
    from ase.io import write
    from ase.data import reference_states,atomic_numbers
    symbol=surface[0].symbol
    crystalstructure=reference_states[atomic_numbers[symbol]]['symmetry']
    if symbol in ['Fe','Co','Ni']:
        spinpol=True
        mixer = MixerFull(beta=0.02,
              nmaxold=5,
              weight=100.0)
        for atom in surface:
            atom.magmom=1.4
    else:
        spinpol=False
        mixer = Mixer(beta=0.05,
              nmaxold=5,
              weight=100.0)
    mol=buildmol(molecule)
    mol.rotate(90,'x')
    adsorbate_indices=ads_map[crystalstructure][molecule]
    mol.positions+=surface.positions[adsorbate_indices[2]]+[0,0,3]
    surface=surface+mol
    surface.calc=GPAW(**calculator,spinpol=spinpol,mixer=mixer)
    qn = BFGS(surface,logfile='%s_%s.log'%(symbol,molecule))
    qn.run(fmax = 0.05)
    write(f'{symbol}_{molecule}_final.traj'),
    return surface


@tb.workflow
class ParametrizableEndStateWorkflow:
    symbol = tb.var()
    calculator = tb.var()

    @tb.task
    def surface(self):
        return tb.node('getsurface',symbol=self.symbol)

    @tb.subworkflow
    def H2(self):
        return EndStateWorkflow(surface=self.surface,molecule='H2',calculator=self.calculator)
    @tb.subworkflow
    def O2(self):
        return EndStateWorkflow(surface=self.surface,molecule='O2',calculator=self.calculator)
    @tb.subworkflow
    def N2(self):
        return EndStateWorkflow(surface=self.surface,molecule='N2',calculator=self.calculator)
    @tb.subworkflow
    def CO(self):
        return EndStateWorkflow(surface=self.surface,molecule='CO',calculator=self.calculator)
    @tb.subworkflow
    def Cl2(self):
        return EndStateWorkflow(surface=self.surface,molecule='Cl2',calculator=self.calculator)
    @tb.subworkflow
    def F2(self):
        return EndStateWorkflow(surface=self.surface,molecule='F2',calculator=self.calculator)


@tb.dynamical_workflow_generator_task
def parametrize_endstate_workflow(calculator):
    material_symbols = [
        'Al', 'Fe','Co','Ni','Cu','Zn','Ru','Rh','Pd','Ag','W','Re','Os','Ir','Pt','Au'
    ]

    for symbol in material_symbols:
        yield f'BEEF-{symbol}', ParametrizableEndStateWorkflow(
            symbol=symbol, calculator=calculator)
