import rasterio as rio
import numpy as np
import fiona
from rasterio.windows import Window
from rasterio import features
from osgeo import gdal
from osgeo import osr

import fiona.crs
from operator import itemgetter
import itertools
from multiprocessing import Pool
from typing import Union
from fct.algorithms.transform import *
from scipy.spatial import KDTree
from shapely.geometry import mapping, Polygon, Point, MultiPoint, multipolygon
from rasterio.features import shapes
from osgeo import ogr
import pandas as pd
import geopandas as gdp
from math import ceil
from collections import OrderedDict
from rasterio import mask



def generate_seeds(feature,type='Poly'):

    if type =='Poly-s':

        for point in feature['geometry']['coordinates']:

            p1 = Polygon(point)
            x, y = p1.centroid.x ,p1.centroid.y

            yield (x, y)
    elif type =='Poly-g':
        for point in feature:

            p1 = point
            x, y = p1.centroid.x ,p1.centroid.y

            yield (x, y)
    else:
        for point in feature['geometry']['coordinates']:

            x, y = point[:2]
            yield (x, y)

def extract_pointxy(name_shp,type_geo='Poly'):

    with fiona.open(name_shp) as ft:
        list_point = [s for feat in ft for s in generate_seeds(feat, type=type_geo)]

    return list_point

def create_grid(xmin,ymin,xmax,ymax,lenght,width,p='General'):

    cols = list(np.arange(np.floor(xmin), np.ceil(xmax), width)) #list(range(int(np.floor(xmin)), int(np.ceil(xmax)), int(width)))
    rows =  list(np.arange(np.floor(ymin), np.ceil(ymax), width)) #list(range(int(np.floor(ymin)), int(np.ceil(ymax)), int(lenght)))

    rows.reverse()
    polygons = []
    for x in cols:
        for y in rows:
            polygons.append(Polygon([(x, y), (x + width, y), (x + width, y - lenght), (x, y - lenght)]))

    if p == 'General':
        centroide_point = [(p1.centroid.x, p1.centroid.y) for p1 in polygons]
    else:
        centroide_point = [(p1.bounds[2], p1.bounds[1]) for p1 in polygons]
    return polygons, centroide_point

def select_pixel_river(river_points, centroide_grid, src_dem, parm):
    kd_tree = KDTree(centroide_grid)

    indexes = kd_tree.query_ball_point(river_points, parm)
    seeds = [(centroide_grid[j][0], centroide_grid[j][1]) for i in range(len(indexes)) for j in indexes[i]]

    meta = src_dem.meta
    # raster_band = raster.read(1)  # read raster band
    pixeles_coord = []
    for i in seeds:
        rowcol = rio.transform.rowcol(meta['transform'], xs=i[0], ys=i[1])
        w = src_dem.read(1, window=Window(rowcol[1], rowcol[0], 1, 1))
        if len(w.tolist()) !=0:
            pixeles_coord.append([i[0], i[1], w.tolist()[0][0]])

    return pixeles_coord


def select_points_buffer(z_value, centroids_point, Meta_raster_dem, src_raster,Threshold):
    select_pixel_by_point = []
    for h in centroids_point:
        xs = h[0]
        ys = h[1]
        rowcol2 = rio.transform.rowcol(Meta_raster_dem, xs=xs, ys=ys)
        w2 = src_raster.read(1, window=Window(rowcol2[1], rowcol2[0], 1, 1))
        row, col = rowcol2
        if len(w2.tolist()) != 0:
            try:
                zs = w2.tolist()[0][0]
            except:
                a = None
        tem_op = zs - z_value     #abs(z_value - zs)
        if  tem_op < Threshold and tem_op!= 0 and tem_op>0 :  # <
            select_pixel_by_point.append([xs, ys, zs, row, col])
    return select_pixel_by_point


def buffer_river_points(xyz_river_points,src_dem,distance_buffer,Threshold):

    lng, wgh = src_dem.res
    meta = src_dem.meta
    pixeles_buffer = []
    contar = 1
    print('Points: ',len(xyz_river_points))
    for j in xyz_river_points:
        print(len(xyz_river_points)-contar)
        z0 = j[2]
        x0 = j[0]
        yo = j[1]
        test = Point(j[0], j[1])


        buf = test.buffer(distance_buffer, cap_style=3)

        xmin, ymin, xmax, ymax = buf.bounds

        p_grids,centroids  = create_grid(xmin, ymin, xmax, ymax, lng, wgh,p='other')

        result = select_points_buffer(z0 , centroids, meta['transform'], src_dem,Threshold)

        pixeles_buffer.extend(result )

        contar += 1
    return pixeles_buffer

def remove_duplicates(points_buffer):

    data = pd.DataFrame(points_buffer, columns=['x', 'y', 'z','row','col'])
    tmp = OrderedDict()
    for point3 in zip(data.x, data.y, data.z,data.row,data.col):
        tmp.setdefault(point3[:2], point3)

    mypoints = tmp.values()


    result = [(pp[0], pp[1], pp[2], pp[3], pp[4]) for pp in mypoints]
    return result

def export_result_raster_and_shp(array_points,cols,rows,raster_reference, output_raster,output_shp,export_Raster=True,export_shp=True):


    zo = pd.DataFrame(np.zeros((cols,rows)))


    for xyz in array_points:
        try:
            zo[xyz[4]].loc[xyz[3]]= xyz[2]
        except:
            print(None)


    zi = np.array(zo)


    if export_Raster == True:
        d = gdal.Open(raster_reference)

        pixelWidth = d.GetGeoTransform()[1]
        pixelHeight = d.GetGeoTransform()[5]
        cols = zi.shape[1]
        rows = zi.shape[0]
        originX = d.GetGeoTransform()[0]
        originY = d.GetGeoTransform()[3]
        driver = gdal.GetDriverByName('GTiff')
        outRaster = driver.Create(output_raster, cols, rows, 1, gdal.GDT_Float32)
        outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
        outband = outRaster.GetRasterBand(1)
        outband.WriteArray(zi)
        outRasterSRS = osr.SpatialReference()
        code_epsg =osr.SpatialReference(wkt= d.GetProjection()).GetAttrValue('AUTHORITY',1)
        outRasterSRS.ImportFromEPSG(int(code_epsg))  # Change EPSG 32618
        outRaster.SetProjection(outRasterSRS.ExportToWkt())
        outband.FlushCache()

    else:
        print(None)

    if export_shp == True:
        mask = None
        with rio.Env():
            with rio.open(output_raster) as src2:
                temp = src2.read(1) # first band
                image = np.where(zi != 0,1, 0)
                results = ({'properties': {'raster_val': v}, 'geometry': s}
                for i, (s, v)
                in enumerate(
                    shapes(image, mask=mask, transform=src2.transform)))

        geoms = list(results)
        gpd_polygonized_raster = gdp.GeoDataFrame.from_features(geoms, crs=src2.meta['crs'])

        gpd_polygonized_raster.to_file(output_shp)
    else:
        print(None)


path_in = r'/'

path_out = r'/'
version = 'v3'
distance =200
threshold = 3
name_drenajes = path_in+'Drenajes_10m_T400.shp'
name_dem =path_in +'dem_10.tif'


outRuta_r = path_out+'Valley'+str(distance/1000)+ 'km_'+str(threshold)+'m_'+version+'.tiff'
outRuta_shp = path_out+'Valley'+str(distance/1000)+ 'km_'+str(threshold)+'m_'+version+'.shp'

points_River = extract_pointxy(name_drenajes,'line')

#tiles = extract_pointxy(path_shp+'grid_arcgis.shp','Poly-s')

# load tif
src = rio.open(name_dem)
pixelH, pixelW = src.res
cols, rows = src.shape
x_min, y_min, x_max, y_max = src.bounds

poly_grids_g,centroids_tiles = create_grid(x_min, y_min, x_max, y_max,pixelH, pixelW)

xyz_river = select_pixel_river(points_River,centroids_tiles,src,pixelH/2)

buffer_points =buffer_river_points(xyz_river,src,distance,threshold)

result_points = remove_duplicates(buffer_points)

export_result_raster_and_shp(array_points=result_points,cols=cols,rows=rows,raster_reference=name_dem,
                             output_raster=outRuta_r,output_shp=outRuta_shp,export_Raster=True,export_shp=True)
