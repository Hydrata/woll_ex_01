#----------------------------------------------------------------------------------------------------------------------------------------------------
# IMPORT NECESSARY MODULES
#----------------------------------------------------------------------------------------------------------------------------------------------------
print 'Setting up Simulation:- Importing Modules'
import anuga, numpy, time, os, glob
from anuga.operators.rate_operators import Polygonal_rate_operator
from anuga import file_function, Polygon_function, read_polygon, create_mesh_from_regions, Domain
import anuga.utilities.quantity_setting_functions as qs
import anuga.utilities.spatialInputUtil as su
import getpass

if __name__ != '__main__':
    from ..database_utils import write_percentage_complete, query_database


def start_woll_ex_01(run_id, Runs, scenario_name, Scenario, session):
    print "start_woll_ex_01 has fired with run_id: %s" % run_id
    base_dir = os.getcwd() + '/tasks/woll_ex_01'
    print "base_dir is: %s" % base_dir

    if run_id == 'local_run':
        base_dir = os.getcwd()

    basename = base_dir + '/inputs/DEM/dem_existing_01'
    outname = run_id
    meshname = base_dir + '/outputs/' + run_id + '.msh'

    # TODO: The database should have already passed the right input files and downloaded them.
    # TODO: A better strategy is to simply import all the files the inputs directories. They should already be correct.
    if run_id == 'local_run':
        bounding_polygon_filename = '%s/inputs/bounding_polygon/%s.shp' % (
            base_dir, 'bdy_02')
    else:
        bounding_polygon_filename = '%s/inputs/bounding_polygon/%s.shp' % (
            base_dir, query_database('bounding_polygon', run_id, Runs, session)[8:])

    print 'bounding_polygon_filename: %s' % bounding_polygon_filename
    bounding_polygon = su.read_polygon(bounding_polygon_filename)
    print 'bounding_polygon: %s' % bounding_polygon

    if run_id == 'local_run':
        interior_holes_filename = '%s/inputs/interior_holes/%s.shp' % (
            base_dir, 'interior_holes_01')
    else:
        interior_holes_filename = '%s/inputs/interior_holes/%s.shp' % (
            base_dir, query_database('interior_holes', run_id, Runs, session)[8:])

    try:
        interior_holes = []
        new_holes = su.read_polygon(interior_holes_filename)
        print 'new_holes: %s' % new_holes
        interior_holes.append(new_holes)
        print 'interior_holes: %s' % interior_holes
    except AssertionError:
        print 'warning: no interior_holes found.'
        interior_holes = None

    # Hack to see how interior_holes should look:
    # interior_holes = [
    #     [
    #         [ 301724.25275325624, 6186056.2254157485 ],
    #         [ 301761.49305447337, 6186006.8056864319 ],
    #         [ 301750.07182667044, 6185997.1975306049 ],
    #         [ 301711.02454357193, 6186045.4673137888 ],
    #         [ 301724.25275325624, 6186056.2254157485 ]
    #     ],
    #     [
    #         [301619.05307519296, 6186067.0245310739],
    #         [301632.49743710633, 6186053.2145431405],
    #         [301623.27550651436, 6186042.9487720868],
    #         [301609.52274617978, 6186055.4096075678],
    #         [301619.05307519296, 6186067.0245310739]
    #     ]
    # ]

    create_mesh_from_regions(bounding_polygon,
        boundary_tags={'south': [0], 'east': [1], 'north': [2], 'west': [3]},
        maximum_triangle_area=50,
        interior_regions=None,
        interior_holes=interior_holes,
        filename=meshname,
        use_cache=False,
        verbose=True)

    domain = Domain(meshname, use_cache=False, verbose=True)
    domain.set_name(outname)
    domain.set_datadir(base_dir + '/outputs')
    print domain.statistics()

    poly_fun_pairs = [['Extent', basename + '.tif']]
    topography_function = qs.composite_quantity_setting_function(
        poly_fun_pairs,
        domain,
        nan_treatment='exception',
    )
    domain.set_quantity('friction', 0.035)
    domain.set_quantity('stage', 0.0)
    domain.set_quantity('elevation', topography_function, verbose=True, alpha=0.99)
    # house_addition_function = qs.composite_quantity_setting_function(house_addition_poly_fun_pairs, domain)
    # domain.add_quantity('elevation', house_addition_function, location='centroids')
    domain.set_minimum_storable_height(0.005)

    #----------------------------------------------------------------------------------------------------------------------------------------------------
    # APPLY RAINFALL
    #----------------------------------------------------------------------------------------------------------------------------------------------------

    Rainfall_Gauge_directory = base_dir + '/inputs/Rainfall/Gauge/'
    for filename in os.listdir(Rainfall_Gauge_directory):
        # only process shapefiles, and only then the .shp:
        if '.shp' in filename:
            Gaugefile = Rainfall_Gauge_directory + filename
            Rainfile = base_dir + '/inputs/Rainfall/rain/rain.tms'
            polygon = su.read_polygon(Gaugefile)
            rainfall = anuga.file_function(Rainfile, quantities='rate')
            op1 = Polygonal_rate_operator(domain, rate=rainfall, factor=1.0e-3, polygon=polygon, default_rate=0.0)

    #----------------------------------------------------------------------------------------------------------------------------------------------------
    # SETUP BOUNDARY CONDITIONS
    #----------------------------------------------------------------------------------------------------------------------------------------------------

    print 'Available boundary tags', domain.get_boundary_tags()

    Br = anuga.Reflective_boundary(domain)
    Bd = anuga.Dirichlet_boundary([0.0,0.0,0.0])

    domain.set_boundary({'interior': Br, 'exterior': Bd, 'west': Bd, 'south': Bd, 'north': Bd, 'east': Bd})

    #----------------------------------------------------------------------------------------------------------------------------------------------------
    # EVOLVE SYSTEM THROUGH TIME
    #----------------------------------------------------------------------------------------------------------------------------------------------------

    t0 = time.time()

    yieldstep = 60.0
    finaltime = 7200.0

    for t in domain.evolve(yieldstep, finaltime):
        domain.write_time()
        percentage_complete = round(domain.time/domain.finaltime, 3)*100
        if run_id != 'local_run':
            write_percentage_complete(run_id, Runs, scenario_name, Scenario, session, percentage_complete)

    print "Done. Nice work."

if __name__ == "__main__":
    # TODO: parse argv for local development
    start_woll_ex_01('local_run', Runs='Runs', session='local_session', Scenario='Scenario', scenario_name='local_scenario')
