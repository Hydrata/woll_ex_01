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
from ..database_utils import write_percentage_complete


def start_woll_ex_01(run_id, Runs, session):
    print "start_woll_ex_01 has fired with run_id: %s" % run_id
    base_dir = os.getcwd() + '/tasks/woll_ex_01'
    print "base_dir is: %s" % base_dir

    #print "linux user is: " + getpass.getuser()
    #----------------------------------------------------------------------------------------------------------------------------------------------------
    #  ADD CATCHMENT INFORMATION HERE
    #----------------------------------------------------------------------------------------------------------------------------------------------------

    # CatchmentDictionary = {
    #     'Model/Bdy/Fine.csv':5
    #     }
    #----------------------------------------------------------------------------------------------------------------------------------------------------
    # FILENAMES, MODEL DOMAIN and VARIABLES
    #----------------------------------------------------------------------------------------------------------------------------------------------------

    basename = base_dir + '/inputs/DEM/dem_existing_01'
    outname = run_id
    meshname = base_dir + '/outputs/' + run_id + '.msh'

    #----------------------------------------------------------------------------------------------------------------------------------------------------
    # ENTER DOMAIN COORDINATES
    #----------------------------------------------------------------------------------------------------------------------------------------------------
    # W=301300
    # N=6186300
    # E=302000
    # S=6185700

    #------------------------------------------------------------------------------
    # CREATING MESH
    #------------------------------------------------------------------------------

    #bounding_polygon = [[W, S], [E, S], [E, N], [W, N]]
    #interior_regions = anuga.read_polygon_dir(CatchmentDictionary, 'Model/Bdy')
    bounding_polygon = su.read_polygon(base_dir + '/inputs/Bdy/bdy_02.shp')
    create_mesh_from_regions(bounding_polygon,
        boundary_tags={'south': [0], 'east': [1], 'north': [2], 'west': [3]},
        maximum_triangle_area=50,
        interior_regions=None,
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
            op1 = Polygonal_rate_operator(domain, rate=rainfall, factor=1.0e-3, polygon=polygon, default_rate = 0.0)

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
    finaltime = 3600.0
    for t in domain.evolve(yieldstep, finaltime):
        domain.write_time()
        percentage_complete = round(domain.time/domain.finaltime, 3)*100
        write_percentage_complete(run_id, Runs, session, percentage_complete)

    print "Done. Nice work."

if __name__ == "__main__":
    # TODO: parse argv for local development
    run_id = 'local_'
    start_woll_ex_01(run_id)