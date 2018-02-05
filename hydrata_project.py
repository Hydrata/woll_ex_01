#----------------------------------------------------------------------------------------------------------------------------------------------------
# IMPORT NECESSARY MODULES
#----------------------------------------------------------------------------------------------------------------------------------------------------
print 'Setting up Simulation:- Importing Modules'
import anuga, numpy, time, os, glob
from anuga.operators.rate_operators import Polygonal_rate_operator
from anuga import file_function, Polygon_function, read_polygon, create_mesh_from_regions, Domain
from anuga import Inlet_operator
import anuga.utilities.quantity_setting_functions as qs
import anuga.utilities.spatialInputUtil as su
import json
from osgeo import ogr

if __name__ != '__main__':
    from ..database_utils import write_percentage_complete, query_database


def start_sim(run_id, Runs, scenario_name, Scenario, session):
    print "start_sim has fired with run_id: %s" % run_id
    base_dir = os.getcwd() + '/base_dir/%s/' % run_id
    project_name = run_id[:10]
    print "base_dir is: %s" % base_dir

    if run_id == 'local_run':
        base_dir = os.getcwd()

    outname = run_id
    meshname = base_dir + 'outputs/' + run_id + '.msh'

    def get_filename(data_type, file_type):
        files = os.listdir('%sinputs/%s' % (base_dir, data_type))
        filename = '%sinputs/%s/%s' % (base_dir,
                                       data_type,
                                       [f for f in files if f[-4:] == file_type][0]
                                       )
        return filename

    bounding_polygon_filename = get_filename('bounding_polygon', '.shp')
    elevation_data_filename = get_filename('elevation_data', '.tif')
    try:
        structures_filename = get_filename('structures', '.shp')
    except OSError as e:
        structures_filename = None
    try:
        rain_data_filename = get_filename('rain_data', '.shp')
    except OSError as e:
        rain_data_filename = None
    try:
        inflow_data_filename = get_filename('inflow_data', '.shp')
    except OSError as e:
        inflow_data_filename = None
    try:
        friction_data_filename = get_filename('friction_data', '.shp')
    except OSError as e:
        friction_data_filename = None

    print 'bounding_polygon_filename: %s' % bounding_polygon_filename
    print 'structures_filename: %s' % structures_filename
    print 'rain_data_filename: %s' % rain_data_filename
    print 'inflow_data_filename: %s' % inflow_data_filename
    print 'friction_data_filename: %s' % friction_data_filename
    print 'elevation_data_filename: %s' % elevation_data_filename
    bounding_polygon = su.read_polygon(bounding_polygon_filename)
    print 'bounding_polygon: %s' % bounding_polygon

    #----------------------------------------------------------------------------------------------------------------------------------------------------
    # SETUP STRUCTURES
    #----------------------------------------------------------------------------------------------------------------------------------------------------

    if structures_filename:
        structures = []
        print 'processing structures from :%s' % structures_filename

        ogr_shapefile = ogr.Open(structures_filename)
        ogr_layer = ogr_shapefile.GetLayer(0)
        ogr_layer_feature = ogr_layer.GetNextFeature()
        while ogr_layer_feature:
            structure = json.loads(ogr_layer_feature.GetGeometryRef().ExportToJson())['coordinates'][0]
            structures.append(structure)
            ogr_layer_feature = None
            ogr_layer_feature = ogr_layer.GetNextFeature()

        print 'structures: %s' % structures
    else:
        print 'warning: no structures found.'
        structures = None

    #----------------------------------------------------------------------------------------------------------------------------------------------------
    # SETUP BOUNDARY CONDITIONS
    #----------------------------------------------------------------------------------------------------------------------------------------------------

    ogr_shapefile = ogr.Open(bounding_polygon_filename)
    ogr_layer = ogr_shapefile.GetLayer(0)
    ogr_layer_definition = ogr_layer.GetLayerDefn()
    print 'ogr_layer_definition.GetGeomType: %s' % ogr_layer_definition.GetGeomType()
    print 'ogr_layer: %s' % ogr_layer
    print 'dir(ogr_layer): %s' % dir(ogr_layer)
    bdy_index = 0
    bdy_tags = {}
    bdy = {}

    ogr_layer_feature = ogr_layer.GetNextFeature()
    while ogr_layer_feature:
        print 'bdy_index: %s' % bdy_index
        boundary_tag_index = int(ogr_layer_feature.GetField('id'))
        boundary_tag_key = ogr_layer_feature.GetField('bdy_tag_k')
        boundary_tag_value = ogr_layer_feature.GetField('bdy_tag_v')
        bdy_tags[boundary_tag_key] = [boundary_tag_index * 2, boundary_tag_index * 2 + 1]
        bdy[boundary_tag_key] = boundary_tag_value
        geom = ogr_layer_feature.GetGeometryRef().GetPoints()
        print 'boundary_tag_key: %s' % boundary_tag_key
        print 'boundary_tag_value: %s' % boundary_tag_value
        print 'geom: %s' % geom
        ogr_layer_feature = None
        ogr_layer_feature = ogr_layer.GetNextFeature()
        bdy_index = bdy_index + 1
        print 'bdy_tags: %s' % bdy_tags

    print 'bdy: %s' % bdy
    print 'create_mesh_from_regions'

    create_mesh_from_regions(
        bounding_polygon,
        boundary_tags=bdy_tags,
        maximum_triangle_area=2,
        interior_regions=None,
        interior_holes=structures,
        filename=meshname,
        use_cache=False,
        verbose=True
    )

    domain = Domain(meshname, use_cache=False, verbose=True)
    domain.set_name(outname)
    domain.set_datadir(base_dir + '/outputs')
    print domain.statistics()
    poly_fun_pairs = [['Extent', elevation_data_filename.encode("utf-8")]]
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

    print 'APPLY RAINFALL'
    if rain_data_filename:
        ogr_shapefile = ogr.Open(rain_data_filename)
        ogr_layer = ogr_shapefile.GetLayer(0)
        ogr_layer_definition = ogr_layer.GetLayerDefn()
        rainfall = 0
        if project_name == 'woll_ex_01':
            rainfile = '/home/ubuntu/anuganode/tasks/woll_ex_01/inputs/Rainfall/rain/rain.tms'
            rainfall = anuga.file_function(rainfile, quantities='rate')
        else:
            ogr_layer_feature = ogr_layer.GetNextFeature()
            while ogr_layer_feature:
                rainfall = ogr_layer_feature.GetField('rate_mm_hr')
                ogr_layer_feature = None
                ogr_layer_feature = ogr_layer.GetNextFeature()
        polygon = su.read_polygon(rain_data_filename)
        Polygonal_rate_operator(domain, rate=rainfall, factor=1.0e-3, polygon=polygon, default_rate=0.0)

    #----------------------------------------------------------------------------------------------------------------------------------------------------
    # APPLY INFLOWS
    #----------------------------------------------------------------------------------------------------------------------------------------------------

    print 'APPLY INFLOWS'
    if inflow_data_filename:
        ogr_shapefile = ogr.Open(inflow_data_filename)
        ogr_layer = ogr_shapefile.GetLayer(0)
        ogr_layer_definition = ogr_layer.GetLayerDefn()
        ogr_layer_feature = ogr_layer.GetNextFeature()
        while ogr_layer_feature:
            in_fixed = ogr_layer_feature.GetField('in_fixed')
            line = ogr_layer_feature.GetGeometryRef().GetPoints()
            Inlet_operator(domain, line, in_fixed, verbose=False)
            ogr_layer_feature = None
            ogr_layer_feature = ogr_layer.GetNextFeature()

    #----------------------------------------------------------------------------------------------------------------------------------------------------
    # APPLY BOUNDARY CONDITIONS
    #----------------------------------------------------------------------------------------------------------------------------------------------------

    print 'APPLY BOUNDARY CONDITIONS'
    print 'Available boundary tags', domain.get_boundary_tags()

    Br = anuga.Reflective_boundary(domain)
    Bd = anuga.Dirichlet_boundary([0.0, 0.0, 0.0])
    Bt = anuga.Transmissive_boundary(domain)

    print 'bdy1: %s' % bdy
    for key, value in bdy.iteritems():
        if value == 'Br':
            bdy[key] = Br
        elif value == 'Bd':
            bdy[key] = Bd
        elif value == 'Bt':
            bdy[key] = Bt
        else:
            print 'No matching boundary condition exists - please check shapefile input'

    print 'bdy2: %s' % bdy
    bdy['interior'] = Br
    print 'bdy3: %s' % bdy

    domain.set_boundary(bdy)
    # domain.set_boundary({'interior': Br, 'exterior': Bd, 'west': Bd, 'south': Bd, 'north': Bd, 'east': Bd})

    #----------------------------------------------------------------------------------------------------------------------------------------------------
    # EVOLVE SYSTEM THROUGH TIME
    #----------------------------------------------------------------------------------------------------------------------------------------------------

    print 'EVOLVE SYSTEM THROUGH TIME'
    t0 = time.time()

    yieldstep = 60.0
    finaltime = 3600.0

    for t in domain.evolve(yieldstep, finaltime):
        domain.write_time()
        percentage_complete = round(domain.time/domain.finaltime, 3)*100
        if run_id != 'local_run':
            write_percentage_complete(run_id, Runs, scenario_name, Scenario, session, percentage_complete)

    print "Done. Nice work."

if __name__ == "__main__":
    # TODO: parse argv for local development
    start_woll_ex_01('local_run', Runs='Runs', session='local_session', Scenario='Scenario', scenario_name='local_scenario')
