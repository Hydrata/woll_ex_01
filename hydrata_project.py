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
import getpass
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

    # TODO: The database should have already passed the right input files and downloaded them.
    # TODO: A better strategy is to simply import all the files the inputs directories. They should already be correct.
    # if run_id == 'local_run':
    #     bounding_polygon_filename = '%s/inputs/bounding_polygon/%s.shp' % (
    #         base_dir, 'bdy_02')
    #     interior_holes_filename = '%s/inputs/interior_holes/%s.shp' % (
    #         base_dir, 'interior_holes_01')
    #     elevation_data_filename = '%s/inputs/elevation_data/%s.shp' % (
    #         base_dir, 'elevation_data')
    # else:
    bounding_polygon_filename = '%sinputs/bounding_polygon/%s.shp' % (
        base_dir, query_database('bounding_polygon', run_id, Runs, session)[8:])
    interior_holes_filename = '%sinputs/interior_holes/%s.shp' % (
        base_dir, query_database('interior_holes', run_id, Runs, session)[8:])
    rain_data_filename = '%sinputs/rain_data/%s.shp' % (
        base_dir, query_database('rain_data', run_id, Runs, session)[8:])
    friction_data_filename = '%sinputs/friction_data/%s.shp' % (
        base_dir, query_database('friction_data', run_id, Runs, session)[8:])
    elevation_data_filename = '%sinputs/elevation_data/%s.tif' % (
        base_dir, query_database('elevation_data', run_id, Runs, session)[8:])

    print 'bounding_polygon_filename: %s' % bounding_polygon_filename
    print 'interior_holes_filename: %s' % interior_holes_filename
    print 'rain_data_filename: %s' % rain_data_filename
    print 'friction_data_filename: %s' % friction_data_filename
    print 'elevation_data_filename: %s' % elevation_data_filename
    bounding_polygon = su.read_polygon(bounding_polygon_filename)
    print 'bounding_polygon: %s' % bounding_polygon

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

    print 'create_mesh_from_regions'
    create_mesh_from_regions(bounding_polygon,
        boundary_tags={'south': [0], 'east': [1], 'north': [2], 'west': [3]},
        maximum_triangle_area=50,
        interior_regions=None,
        interior_holes=interior_holes,
        filename=meshname,
        use_cache=False,
        verbose=True)

    print 'domain'
    domain = Domain(meshname, use_cache=False, verbose=True)
    domain.set_name(outname)
    domain.set_datadir(base_dir + '/outputs')
    print domain.statistics()

    print 'create poly_fun_pairs'
    poly_fun_pairs = [['Extent', elevation_data_filename.encode("utf-8")]]
    print 'set poly_fun_pairs: %s' % poly_fun_pairs
    print 'create topography_function'
    topography_function = qs.composite_quantity_setting_function(
        poly_fun_pairs,
        domain,
        nan_treatment='exception',
    )
    print 'domain.set_quantity'
    domain.set_quantity('friction', 0.035)
    domain.set_quantity('stage', 0.0)
    print 'domain.set_quantity elevation'
    domain.set_quantity('elevation', topography_function, verbose=True, alpha=0.99)
    # house_addition_function = qs.composite_quantity_setting_function(house_addition_poly_fun_pairs, domain)
    # domain.add_quantity('elevation', house_addition_function, location='centroids')
    domain.set_minimum_storable_height(0.005)

    #----------------------------------------------------------------------------------------------------------------------------------------------------
    # APPLY RAINFALL
    #----------------------------------------------------------------------------------------------------------------------------------------------------

    ogr_shapefile = ogr.Open(rain_data_filename)
    ogr_layer = ogr_shapefile.GetLayer(0)
    ogr_layer_definition = ogr_layer.GetLayerDefn()
    print 'ogr_layer_definition.GetGeomType: %s' % ogr_layer_definition.GetGeomType()

    print 'APPLY RAINFALL'
    if ogr_layer_definition.GetGeomType() == 3:
        # Rainfall_Gauge_directory = base_dir + '/inputs/rainfall_data/'
        # for filename in os.listdir(Rainfall_Gauge_directory):
        #     # only process shapefiles, and only then the .shp:
        #     if '.shp' in filename:
        #        Gaugefile = Rainfall_Gauge_directory + filename
        #
        if project_name == 'woll_ex_01':
            rainfile = '/home/ubuntu/anuganode/tasks/woll_ex_01/inputs/Rainfall/rain/rain.tms'
            rainfall = anuga.file_function(rainfile, quantities='rate')
        else:
            ogr_layer_feature = ogr_layer.GetNextFeature()
            while ogr_layer_feature:
                rain_fixed = ogr_layer_feature.GetField('rain_fixed')
                geom = ogr_layer_feature.GetGeometryRef().GetPoints()
                print 'in_fixed: %s' % rain_fixed
                print 'geom: %s' % geom
                Inlet_operator(domain, geom, rain_fixed, verbose=False)
                ogr_layer_feature = None
                ogr_layer_feature = ogr_layer.GetNextFeature()

        # print 'dir rainfall: %s' % dir(rainfall)
        # print 'rainfall.vertex_coordinates: %s' % rainfall.vertex_coordinates
        # print 'rainfall.quantity_names: %s' % rainfall.quantity_names
        # print 'rainfall.statistics: %s' % rainfall.statistics()
        # print 'rainfall.quantities_range: %s' % rainfall.quantities_range
        # print 'rainfall.precomputed_values: %s' % rainfall.precomputed_values
        # print 'rainfall.index: %s' % rainfall.index
        # print 'rainfall(): %s' % rainfall()
        polygon = su.read_polygon(rain_data_filename)
        Polygonal_rate_operator(domain, rate=rainfall, factor=1.0e-3, polygon=polygon, default_rate=0.0)

    #----------------------------------------------------------------------------------------------------------------------------------------------------
    # APPLY INFLOWS
    #----------------------------------------------------------------------------------------------------------------------------------------------------

    print 'APPLY INFLOWS'
    if ogr_layer_definition.GetGeomType() == 2:
        # print 'ogr_layer_definition.GetGeomFieldDefn: %s' % ogr_layer_definition.GetGeomFieldDefn()
        # print 'ogr_layer_definition.GetGeomFieldCount: %s' % ogr_layer_definition.GetGeomFieldCount()
        # print 'ogr_layer_definition.GetGeomFieldIndex: %s' % ogr_layer_definition.GetGeomFieldIndex()

        print "Name  -  Type  Width  Precision"
        for i in range(ogr_layer_definition.GetFieldCount()):
            fieldName = ogr_layer_definition.GetFieldDefn(i).GetName()
            fieldTypeCode = ogr_layer_definition.GetFieldDefn(i).GetType()
            fieldType = ogr_layer_definition.GetFieldDefn(i).GetFieldTypeName(fieldTypeCode)
            fieldWidth = ogr_layer_definition.GetFieldDefn(i).GetWidth()
            GetPrecision = ogr_layer_definition.GetFieldDefn(i).GetPrecision()
            print fieldName + " - " + fieldType + " " + str(fieldWidth) + " " + str(GetPrecision)

        print 'ogr_layer_definition: %s' % dir(ogr_layer_definition)

        ogr_layer_feature = ogr_layer.GetNextFeature()
        while ogr_layer_feature:
            in_fixed = ogr_layer_feature.GetField('in_fixed')
            line = ogr_layer_feature.GetGeometryRef().GetPoints()
            print 'in_fixed: %s' % in_fixed
            print 'geom: %s' % line
            Inlet_operator(domain, line, in_fixed, verbose=False)
            ogr_layer_feature = None
            ogr_layer_feature = ogr_layer.GetNextFeature()

    #----------------------------------------------------------------------------------------------------------------------------------------------------
    # SETUP BOUNDARY CONDITIONS
    #----------------------------------------------------------------------------------------------------------------------------------------------------

    print 'SETUP BOUNDARY CONDITIONS'
    print 'Available boundary tags', domain.get_boundary_tags()

    Br = anuga.Reflective_boundary(domain)
    Bd = anuga.Dirichlet_boundary([0.0, 0.0, 0.0])

    domain.set_boundary({'interior': Br, 'exterior': Bd, 'west': Bd, 'south': Bd, 'north': Bd, 'east': Bd})

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
