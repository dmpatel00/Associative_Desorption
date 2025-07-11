#
#
#
#
from ase.db import connect
from ase.data import reference_states,atomic_numbers
from os.path import isfile
elements = ['Al','Ni','Cu','Rh','Pd','Ag','Ir','Pt','Au','Co','Zn','Ru','Re','Os', 'Fe','W']
adsorbates=['H','N','O','C','Cl','F']
sitedict={
        'fcc':['ontop','bridge','fcc','hcp'],
        'hcp':['ontop','bridge','fcc','hcp'],
        'bcc':['ontop','longbridge','shortbridge','hollow']
        }
with connect('adsorbates.db') as masterdb:
    for element in elements:
        surfdb=connect(f'tree/surfaces/BEEF-{element}/min_surf/relax/{element}_clean.db')
        try:
            surfrow=surfdb.get(element=element,adsorbate='clean')
        except AssertionError:
                print(f'{element}_surface has mulitple entries')
                energies=[]
                ids=[]
                for row in surfdb.select(element=element,adsorbate='clean'):
                    energies.append(row.energy)
                    ids.append(row.id)
                    surfrow=surfdb.get(id=ids[energies.index(min(energies))])
        crystalstructure=reference_states[atomic_numbers[element]]['symmetry']
        sites=sitedict[crystalstructure]
        masterdb.write(surfrow,element=element,adsorbate='clean',site='None',crystal=crystalstructure)
        for adsorbate in adsorbates:
            if  isfile(f'tree/surfaces/BEEF-{element}/add_{adsorbate}/relax/{element}_{adsorbate}.db'):
                subdb=connect(f'tree/surfaces/BEEF-{element}/add_{adsorbate}/relax/{element}_{adsorbate}.db')
                for site in sites:
                    try:
                        row=subdb.get(element=element,adsorbate=adsorbate,site=site)
                        masterdb.write(row,element=element,adsorbate=adsorbate,site=site,crystal=crystalstructure)
                    except KeyError:
                        print(f'{element}_{adsorbate}_{site} not found')
                    except AssertionError:
                        print(f'{element}_{adsorbate}_{site} has mulitple entries')
                        energies=[]
                        ids=[]
                        for row in subdb.select(element=element,adsorbate=adsorbate,site=site):
                            energies.append(row.energy)
                            ids.append(row.id)
                        row=subdb.get(id=ids[energies.index(min(energies))])
                        masterdb.write(row,element=element,adsorbate=adsorbate,site=site,crystal=crystalstructure)
