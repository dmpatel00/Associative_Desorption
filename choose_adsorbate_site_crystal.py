#
#
#
#
from ase.db import connect
import numpy as np
db=connect('adsorbates.db')
adsorbates=['H','N','O','C','Cl','F']
crystals=['hcp']#,'hcp']#,'bcc']
for adsorbate in adsorbates:
    for crystal in crystals:
        selection=db.select(crystal=crystal,adsorbate=adsorbate)
        data=np.empty((0,2))
        energies=np.array([])
        sites=[]
        for row in selection:
            data=np.append(data,[[row.element, row.site]],axis=0)
            energies=np.append(energies,row.energy)
        elements=np.unique(data[:,[0]])
        for element in elements:
            site_energies=energies[data[:,0]==element]
            min_energy=np.min(site_energies)
            #print(data[energies==min_energy][0])
            sites.append(data[energies==min_energy][0][1])
        print(elements)
        print(crystal, adsorbate, sites) 
        print(max(set(sites),key=sites.count))
