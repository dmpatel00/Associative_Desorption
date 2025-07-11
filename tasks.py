from pathlib import Path
from gpaw import GPAW
import taskblaster as tb
import numpy as np

RPBE_latticeconstants={
        'Ag': 2.1032588528146277,
        'Al': 2.034719160256986,
        'Au': 2.1039660814681618,
        'Cu': 1.8374052933055884,
        'Fe': 1.4342135039949488,
        'Ir': 1.9417086002991968,
        'Ni': 1.7814747787301133,
        'Pd': 1.9920604957429109,
        'Pt': 1.99649547785747,
        'Rh': 1.933116214751228,
        'W': 1.5988082893039466,
        'Co':[ 2.490610126722544, 4.039138843362412],
        'Os':[ 2.7589329992487936, 4.357713272155947],
        'Re':[ 2.787459301514976, 4.4905075709326345],
        'Ru':[ 2.731233112069678, 4.311485469812148],
        'Zn':[ 2.7639362406148056, 4.68600399952026],
        }

def optimize_cell(atoms,symbol,calculator):
    from gpaw import GPAW, Mixer, MixerFull
    from ase.build import bulk
    from ase.data import reference_states,atomic_numbers
    crystalstructure=reference_states[atomic_numbers[symbol]]['symmetry']
    lat_const=RPBE_latticeconstants[symbol]
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
    from ase.io import Trajectory
    traj=Trajectory(f'{symbol}.traj','w')
    eps = 0.05
    energies=np.array([])
    if crystalstructure=='hcp':
        a_consts=np.array([])
        c_consts=np.array([])
        a0=lat_const[0]
        c0=lat_const[1]
        for a in a0 * np.linspace(1 + eps, 1 - eps, 11):
            for c in c0 * np.linspace(1 + eps, 1 - eps, 11):
                a_consts=np.append(a_consts,a)
                c_consts=np.append(c_consts,c)
                atoms=bulk(symbol,crystalstructure,a=a,c=c)
                atoms.calc = GPAW(**calculator,spinpol=spinpol,mixer=mixer)
                energy=atoms.get_potential_energy()
                energies=np.append(energies,energy)
                traj.write(atoms)
        min_index=np.argmin(energies)
        lattice_outa=a_consts[min_index]
        lattice_outc=c_consts[min_index]
        min_atoms=bulk(symbol,crystalstructure,a=lattice_outa,c=lattice_outc)
    else:
        a_consts=np.array([])
        a0=2*lat_const
        for a in a0 * np.linspace(1 + eps, 1 - eps, 11):
            a_consts=np.append(a_consts,a)
            atoms=bulk(symbol,crystalstructure,a=a)
            atoms.calc = GPAW(**calculator,spinpol=spinpol,mixer=mixer)
            energy=atoms.get_potential_energy()
            energies=np.append(energies,energy)
            traj.write(atoms)
        min_index=np.argmin(energies)
        lattice_out=a_consts[min_index]
        min_atoms=bulk(symbol,crystalstructure,a=lattice_out)
    return min_atoms

def asebulk(symbol):
    from ase.build import bulk
    return bulk(symbol)

@tb.workflow
class MaterialsWorkflow:
    atoms = tb.var()
    symbol = tb.var()
    calculator = tb.var()

    @tb.task
    def relax(self):
        return tb.node(
            'optimize_cell',
            atoms=self.atoms,
            symbol=self.symbol,
            calculator=self.calculator)



@tb.workflow
class ParametrizableMaterialsWorkflow:
    symbol = tb.var()
    calculator = tb.var()

    @tb.task
    def atoms(self):
        return tb.node('asebulk', symbol=self.symbol)

    @tb.subworkflow
    def compute(self):
        return MaterialsWorkflow(atoms=self.atoms,symbol=self.symbol, calculator=self.calculator)

   # @tb.task
   #     return tb.node('buildsurface',dimensions=(3,4,3),symbol=self.symbol)
   #     return tb.node('buildsurface',dimensions=(3,4,3),symbol=self.symbol)


@tb.dynamical_workflow_generator_task
def parametrize_materials_workflow(calculator):
    from ase.build import bulk

    material_symbols = [
        'Al', 'Fe','Co','Ni','Cu','Zn','Ru','Rh','Pd','Ag','W','Re','Os','Ir','Pt','Au'
    ]

    for symbol in material_symbols:
        yield f'BEEF-{symbol}', ParametrizableMaterialsWorkflow(
            symbol=symbol, calculator=calculator)

distances={
        'H':3.5,
        'C':3.5,
        'O':3.5,
        'Cl':3.8,
        'N':4.2,
        'F':3.8,
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
        surface=bcc110(symbol,dimensions,a=2*lat_const,orthogonal=True)
    from ase.constraints import FixAtoms
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
        maxiter=1000
        spinpol=True
        mixer = MixerFull(beta=0.05,
              nmaxold=5,
              weight=100.0)
        for atom in surface:
            atom.magmom=2.8
    else:
        maxiter=333
        spinpol=False
        mixer = Mixer(beta=0.02,
              nmaxold=5,
              weight=50.0)
    site='None'
    try:
        surface.calc=GPAW(**calculator,spinpol=spinpol,mixer=mixer,maxiter=maxiter)
        qn = BFGS(surface,logfile='%s_%s_%s.log'%(symbol,adsorbate,site))
        qn.run(fmax = 0.05)
    except RuntimeError:
        surface.calc=GPAW(**calculator,spinpol=spinpol,mixer=mixer,symmetry='off')
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
    from ase.io import write,read
    from gpaw import GPAW, Mixer, MixerFull
    from ase.optimize import BFGS
    from ase.db import connect
    from os.path import isfile
    from gpaw import Davidson
    symbol=surface[0].symbol
    db=connect(f'{symbol}_{adsorbate}.db')
    if symbol in ['Fe','Co','Ni']:
        maxiter=1000
        spinpol=True
        mixer = MixerFull(beta=0.02,
              nmaxold=3,
              weight=100.0)
        for atom in surface:
            atom.magmom=2.8
    else:
        maxiter=500
        spinpol=False
        mixer = Mixer(beta=0.05,
              nmaxold=5,
              weight=100.0)
    sites=surface.info['adsorbate_info']['sites'].keys()
    min_atoms_vec=[]
    for site in sites:
        try:
            db.get(site=site)
            continue
        except:
            pass
        previous_exists=isfile('Traj_%s_%s_%s.traj'%(symbol,adsorbate,site))
        if previous_exists:
            try:
                surf_ads=read('Traj_%s_%s_%s.traj'%(symbol,adsorbate,site))
            except:
                surf_ads=surface.copy()
                add_adsorbate(surf_ads,adsorbate,distances[adsorbate],site)
        else:
            surf_ads=surface.copy()
            add_adsorbate(surf_ads,adsorbate,distances[adsorbate],site)
        surf_ads.calc=GPAW(**calculator,spinpol=spinpol,mixer=mixer,maxiter=maxiter,eigensolver=Davidson(3))
        qn = BFGS(surf_ads,trajectory='Traj_%s_%s_%s.traj'%(symbol,adsorbate,site),logfile='%s_%s_%s.log'%(symbol,adsorbate,site))
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
        return tb.node('buildsurf',dimensions=(3,4,3),symbol=self.symbol)

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

def getsurface(symbol):
    from ase.db import connect
    db=connect('/home/cat/dmapa/gpaw/dissads/adsorbates.db')
    surface=db.get(element=symbol,adsorbate='clean').toatoms()
    return surface

ads_map={
        'fcc':{
            'H2':[16,20,32],
            'N2':[16,20,32],
            'O2':[16,20,32],
            'CO':[16,20,32],
            'CN':[16,20,32],
            'F2':[16,20,32],
            'Cl2':[16,20,32],
            },
        'hcp':{
            'H2':[13,16,31],
            },
        'bcc':{
            'H2':[13,16,31],
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
    def initial_andreas(self):
        return tb.node(
            'initial_catlearn',
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

    @tb.task
    def doNEB(self):
        return tb.node(
            'do_NEB',
            start=self.initial,
            end=self.final,
            molecule=self.molecule,
            calculator=self.calculator)

def do_NEB(start,end,molecule,calculator):
    from gpaw import GPAW, Mixer, MixerFull
    from ase.io import read,write
    import sys,glob
    sys.path.insert(0,'/home/cat/dmapa/venv0/CatLearn')
    from catlearn.optimize.mlneb import MLNEB
    from catlearn.optimize.neb import ImprovedTangentNEB
    initial=read(glob.glob('../initial/*_initial.traj')[0])
    final=read(glob.glob('../final/*_final.traj')[0])
    symbol=initial[0].symbol
    if symbol in ['Fe','Co','Ni']:
        spinpol=True
        mixer = MixerFull(beta=0.02,
              nmaxold=5,
              weight=100.0)
        maxiter=1000
    else:
        spinpol=False
        mixer = Mixer(beta=0.05,
              nmaxold=5,
              weight=100.0)
        maxiter=333
    if molecule in ['O2','CN']:
        spinpol=True
        mixer = MixerFull(beta=0.02,
              nmaxold=5,
              weight=100.0)
        maxiter=500
    n_images=15
    unc_convergence=0.05
    max_unc=0.50
    calc=GPAW(**calculator,spinpol=spinpol,mixer=mixer,maxiter=maxiter)
    mlneb=MLNEB(start=initial,end=final,
                ase_calc=calc,
                interpolation='linear', #changed for Pd
                interpolation_kwargs=dict(mic=False),
                neb_method=ImprovedTangentNEB,
                n_images=n_images,
                full_output=True)
    mlneb.run(fmax=0.05,unc_convergence=unc_convergence,max_unc=max_unc,steps=100,ml_steps=500)

    return initial

def initial_catlearn(surface,molecule,calculator):
    import sys
    sys.path.insert(0,'/home/cat/alyvi/pkg_versions/CatLearn_versions/CatLearn601/')
    from gpaw import GPAW, Mixer, MixerFull
    from ase.optimize import BFGS
    from ase.io import write
    from ase.data import reference_states,atomic_numbers
    from ase.build import molecule as buildmol
    import numpy as np
    from catlearn.activelearning.adsorption import AdsorptionAL
    from ase import Atoms
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
    if molecule in ['O2','CN']:
        spinpol=True
        mixer = MixerFull(beta=0.02,
              nmaxold=5,
              weight=100.0)
    adsorbate_indices=ads_map[crystalstructure][molecule]
    mol=buildmol(molecule)
    ads=Atoms(mol[0].symbol)
    ads2=Atoms(mol[1].symbol)
    calc=GPAW(**calculator,spinpol=spinpol,mixer=mixer)
    bounds = np.array([
        [0.0, 1.0],
        [0.0, 1.0],
        [0.55, 0.8],
        [0.0, 2 * np.pi],
        [0.0, 2 * np.pi],
        [0.0, 2 * np.pi],
        [0.0, 1.0],
        [0.0, 1.0],
        [0.55, 0.8],
        [0.0, 2 * np.pi],
        [0.0, 2 * np.pi],
        [0.0, 2 * np.pi],
        ])
    dyn=AdsorptionAL(
            slab=surface,
            adsorbate=ads,
            adsorbate2=ads2,
            ase_calc=calc,
            unc_convergence=0.02,
            bounds=bounds,
            opt_kwargs={},
            parallel_run=False,
            min_data=3,
            restart=False,
            verbose=True)
    dyn.run(
            fmax=0.05,
            max_unc=0.30,
            steps=100,
            ml_steps=4000)
    write(f'{symbol}_{molecule}_initial.traj',dyn.structures)
    return surface

def initial_state(surface,molecule,calculator):
    from gpaw import GPAW, Mixer, MixerFull
    from ase.optimize import BFGS
    from ase.io import write
    from ase.data import reference_states,atomic_numbers
    from ase.build import molecule as buildmol
    symbol=surface[0].symbol
    crystalstructure=reference_states[atomic_numbers[symbol]]['symmetry']
    if symbol in ['Fe','Co','Ni']:
        spinpol=True
        mixer = MixerFull(beta=0.02,
              nmaxold=5,
              weight=100.0)
        for atom in surface:
            atom.magmom=1.4
        maxiter=1000
    else:
        spinpol=False
        mixer = Mixer(beta=0.05,
              nmaxold=5,
              weight=100.0)
        maxiter=333
    if molecule in ['O2','CN']:
        spinpol=True
        mixer = MixerFull(beta=0.02,
              nmaxold=5,
              weight=100.0)
    adsorbate_indices=ads_map[crystalstructure][molecule]
    mol=buildmol(molecule)
    surface.append(mol[0].symbol)
    surface.positions[-1]=surface.positions[adsorbate_indices[0]]+(0,0,distances[mol[0].symbol])
    surface.append(mol[1].symbol)
    surface.positions[-1]=surface.positions[adsorbate_indices[1]]+(0,0,distances[mol[1].symbol])
    if crystalstructure in ['hcp','bcc']:
        surface.positions[-1]+=(0,1.4,0)
        surface.positions[-2]+=(0,1.4,0)
    surface.calc=GPAW(**calculator,spinpol=spinpol,mixer=mixer,maxiter=maxiter)
    qn = BFGS(surface,logfile='%s_%s.log'%(symbol,molecule))
    qn.run(fmax = 0.05)
    write(f'{symbol}_{molecule}_initial.traj',surface)
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
    if molecule in ['O2','CN']:
        spinpol=True
        mixer = MixerFull(beta=0.02,
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
    write(f'{symbol}_{molecule}_final.traj',surface)
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
    def CN(self):
        return EndStateWorkflow(surface=self.surface,molecule='CN',calculator=self.calculator)
    @tb.subworkflow
    def Cl2(self):
        return EndStateWorkflow(surface=self.surface,molecule='Cl2',calculator=self.calculator)
    @tb.subworkflow
    def F2(self):
        return EndStateWorkflow(surface=self.surface,molecule='F2',calculator=self.calculator)


@tb.dynamical_workflow_generator_task
def parametrize_endstate_workflow(calculator):
    #material_symbols = [
    #    'Al', 'Fe','Co','Ni','Cu','Zn','Ru','Rh','Pd','Ag','W','Re','Os','Ir','Pt','Au'
    #]
    material_symbols = [
        'Al','Ni','Cu','Rh','Pd','Ag','Ir','Pt','Au','Co','Zn','Ru','Re','Os', 'Fe','W'
    ]

    for symbol in material_symbols:
        yield f'BEEF-{symbol}', ParametrizableEndStateWorkflow(
            symbol=symbol, calculator=calculator)
