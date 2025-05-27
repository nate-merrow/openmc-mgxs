import openmc
import os

os.makedirs("output", exist_ok=True)
os.chdir("output")

##############################################
                # Materials #
##############################################

uo2 = openmc.Material(name='uo2 fuel')
uo2.add_nuclide('U235',0.04)
uo2.add_nuclide('U238',0.96)
uo2.add_element('O',2.00)
uo2.set_density('g/cm3',10.4)

zircaloy = openmc.Material(name='zircaloy cladding')
zircaloy.add_element('Zr',0.982)
zircaloy.add_element('Sn',0.015)
zircaloy.add_element('Fe',0.002)
zircaloy.add_element('Cr',0.001)
zircaloy.set_density('g/cm3',6.55)

water = openmc.Material(name='water')
water.add_nuclide('H1',2.00)
water.add_nuclide('O16',1.00)
water.set_density('g/cm3',1.00)
water.add_s_alpha_beta('c_H_in_H2O')

materials = openmc.Materials([uo2, zircaloy, water])
materials.export_to_xml()

##############################################
                # Geometry #
##############################################

fuel_outer_radius = openmc.ZCylinder(r=0.4096)
clad_inner_radius = openmc.ZCylinder(r=0.4179)
clad_outer_radius = openmc.ZCylinder(r=0.4750)

fuel_region = -fuel_outer_radius
gap_region = +fuel_outer_radius & -clad_inner_radius
clad_region = +clad_inner_radius & -clad_outer_radius

fuel = openmc.Cell(1,name='fuel')
fuel.fill = uo2
fuel.region = fuel_region

gap = openmc.Cell(name='air gap')
gap.region = gap_region

clad = openmc.Cell(name='clad')
clad.fill = zircaloy
clad.region = clad_region

pitch = 1.26
box = openmc.model.RectangularPrism(width=pitch, height=pitch, boundary_type='reflective')
type(box)

water_region = +clad_outer_radius & -box

moderator = openmc.Cell(name='water moderator')
moderator.fill = water
moderator.region = water_region

root_universe = openmc.Universe(cells=[fuel, gap, clad, moderator])

geometry = openmc.Geometry()
geometry.root_universe = root_universe
geometry.export_to_xml()

##############################################
        # Multi-Group Cross Sections #
##############################################

import openmc.mgxs

energy_groups = openmc.mgxs.EnergyGroups(group_edges=[1.0e-3, 1.0, 0.1e6, 10e6]) # based on E.E.Lewis page 40

mgxs_library = openmc.mgxs.Library(geometry) # mgxs library creation
mgxs_library.energy_groups = energy_groups
mgxs_library.domain_type = 'cell'
mgxs_library.domains = [fuel] # only record cross sections for the fuel
mgxs_library.mgxs_types = ['total', 'fission', 'absorption', 'scatter matrix']
mgxs_library.build_library()

##############################################
                # Tallies #
##############################################

cell_filter = openmc.CellFilter(fuel) # tally only from fuel

tally = openmc.Tally(1)
tally.filters = [cell_filter]
tally.scores = ['total', 'fission', 'absorption', '(n,gamma)']

rr_tally = openmc.Tally(name='reaction_rates')
rr_tally.scores = ['fission', 'absorption']

tallies = openmc.Tallies([tally, rr_tally])
mgxs_library.add_to_tallies_file(tallies)
tallies.export_to_xml()

##############################################
                # Settings #
##############################################

source = openmc.IndependentSource(space=openmc.stats.Point((0,0,0)))

settings = openmc.Settings()
settings.batches = 100
settings.inactive = 10
settings.particles = 10000
settings.source = source

settings.export_to_xml()

##############################################
                 # Run #
##############################################

openmc.run()

##############################################
             # Post-Processing #
##############################################

import matplotlib.pyplot as plt 
import numpy as np

sp = openmc.StatePoint('statepoint.100.h5') # access the statepoint file in output
mgxs_library.load_from_statepoint(sp) # load the library from the statepoint

for domain in mgxs_library.domains:
    for xs_type in mgxs_library.mgxs_types:
        mgxs = mgxs_library.get_mgxs(domain, xs_type) # create mgxs object from the mgxs library that can be worked with more easily

        print(f"\nCross section type: {xs_type.upper()} for domain ID: {domain.id}")

        if xs_type == 'scatter matrix':
            xs_matrix = mgxs.get_xs()
            
            fig, ax = plt.subplots() # heat map plotting for the scatter matrix
            c = ax.imshow(xs_matrix, cmap='viridis', origin='lower') # create the color map
            fig.colorbar(c, ax=ax)
            ax.set_title(f"Scatter Matrix - Domain {domain.id}")
            ax.set_xlabel("To Group")
            ax.set_ylabel("From Group")
            plt.tight_layout()
            plt.savefig('scattering.png') # same the plot to this png file in output
                
        else:     # for all scalar mgxs types
            try:
                xs_data = mgxs.get_xs(nuclide='total')
            except TypeError:
                xs_data = mgxs.get_xs()
                
            print(f"XS per energy group: {xs_data.shape}")
            print(xs_data)
            
            energy_groups_again = mgxs.energy_groups.group_edges
            group_midpoints = 0.5 * (energy_groups_again[:-1] + energy_groups_again[1:]) # binning for the bar graph
            group_labels = [f"{e/1e6:.2f} MeV" for e in group_midpoints]

            plt.figure()
            plt.bar(range(len(xs_data)), xs_data, tick_label=group_labels)
            plt.xticks(rotation=45, ha='right')
            plt.ylabel("Cross Sections [barns]")
            plt.title(f"{xs_type.title()} XS - Domain {domain.id}")
            plt.tight_layout
            plt.savefig(f"{xs_type.title()}_XS.png") # saving each plot to a png file based on the cross section types name







    



