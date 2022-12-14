import fct.terrain_analysis as ta
import rasterio as rio
import numpy as np
import fiona
from osgeo import gdal
from osgeo import osr
import fct.speedup as ss
import fiona.crs



def array2raster(newRasterfn, array, dir_DEM):
    pixelWidth = gdal.Open(dir_DEM).GetGeoTransform()[1]
    pixelHeight = gdal.Open(dir_DEM).GetGeoTransform()[5]
    cols = array.shape[1]
    rows = array.shape[0]
    originX = gdal.Open(dir_DEM).GetGeoTransform()[0]
    originY = gdal.Open(dir_DEM).GetGeoTransform()[3]
    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(newRasterfn, cols, rows, 1, gdal.GDT_Int32)
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(array)
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromEPSG(32632)  # Change EPSG 32618
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()

def LoadRaster(file_name, showIndicator=False):
    '''Esta funcion carga un raster usando GDAL;
    file_name - contiene la ruta completa hacia el archivo raster;
    showIndicator - Activa o desactiva el mensaje de carga de archivo
    grid - es la grilla raster que regresa la funcion.
    la grilla que se regresa es un arreglo de numpy
    En el codigo esta funcion sirve para cargar los grd de surfer
    de precipitacion y el DEM'''

    g = gdal.Open(file_name, gdal.GA_ReadOnly)
    g.GetRasterBand(1) # We read the first band allways
    cols = g.RasterXSize
    nodatavalue = g.GetRasterBand(1).GetNoDataValue()
    rows = g.RasterYSize
    Cell_size = g.GetGeoTransform()[1]
    grid = g.ReadAsArray(0, 0, cols, rows)
    if showIndicator == True:
        print ('Cargando el archivo==> ' + file_name)
    else:
        pass
    return grid, nodatavalue,Cell_size

def fill_dem(array, nodata):
    data = np.float32(array)  # trasformacion de datos a valores numericos float 32 bits
    try:
        fill_dem = ta.fillsinks(data, nodata)  # proceso de llenado de celdas vacias
    except:
        nodata = -32768
        fill_dem = ta.fillsinks(data, -32768)
    print('Fill ','Done')
    return fill_dem, nodata

def flow_direction(array, nodata):
    flow = ta.flowdir(array, nodata, flow=None)  # proceso de estimacion de flujo

    mask, labels = ta.resolve_flat(array, flow)  # mascara de flujo

    flowdir = ta.flat_mask_flowdir(mask, flow, labels)  # Estimacion de direccion de flujo
    print('Flow Direction ', 'Done')
    return flowdir

def flow_acumulation(array_flow):
    print('Flow Acumulation ', 'Done')
    return ta.flow_accumulation(array_flow)

def create_threshold(array_dir, array_acu, size_cell, threshold):
    CellA = (float(size_cell) * float(size_cell)) / 500000  # 1000000
    nbCell = float(threshold) / CellA  # estimacion del limite numero de celdas

    # ThresholdStep = np.where(flow_acu>=nbCell,flow_acu,np.NaN)
    ThresholdStep = np.int16(array_acu > nbCell)  # exportar
    feature_obj = ss.stream_to_feature(ThresholdStep,
                                       np.int16(array_dir))  # generacion de objeto que contiene drenajes fragmentados
    print('Objet Threshold ', 'Done')
    return feature_obj

def export_shp(array_flowac,obj_features,name_out_shp,name_dem):
    ref_raster = rio.open(name_dem)
    epsg_out = int(ref_raster.crs['init'][5:])
    # salidas shapefile y raster
    schema = {'geometry': 'LineString',
        'properties': [('GID', 'int'),
            ('HEAD', 'int:1'),
            ('ROW', 'int:4'),
            ('COL', 'int:4')]}

    driver = 'ESRI Shapefile'
    crs = fiona.crs.from_epsg(epsg_out) #sistema de coordenada de salida
    options = dict(driver=driver, crs=crs, schema=schema)
    row, col = array_flowac.shape[0],array_flowac.shape[1] #configuracion de shp de salida
    with fiona.open(name_out_shp, 'w', **options) as dst:
        for current, (segment, head) in enumerate (obj_features ):
            coords = ta.pixeltoworld(np.fliplr(np.int32(segment)), ref_raster.transform, gdal=False)
            dst.write({
                'type': 'Feature',
                'geometry': {'type': 'LineString', 'coordinates': coords},
                'properties': {
                    'GID': current,
                    'HEAD': 1 if head else 0,
                    'ROW': row,
                    'COL': col
                }
            })
    print('Export Shapefile ', 'Done')




path_in = r''
path_out = r''


name_dem ='dem_10.tif' # 'dem.tif'#'Dem_EX_1.tif' # nombre raster
Suflix = '10m'

Threshold = 400

data, X_No_Data, cellZ = LoadRaster(path_in+name_dem ) # leer dem de entrada

array_fill,X_No_Data = fill_dem(data,X_No_Data)

flowd_array = flow_direction(array_fill,X_No_Data)

flow_acu =flow_acumulation(flowd_array)

feature_ob = create_threshold(flowd_array, flow_acu, cellZ, Threshold)

export_shp(flow_acu,feature_ob,path_out+'Drenajes_'+Suflix+'_T'+str(Threshold),path_in+name_dem)

