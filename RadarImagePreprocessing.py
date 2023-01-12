import os
import datetime


# snappy operations

def read_product(src, type='.dim'):
    """
    Wczytaj produkt. Domyślny format to '.dim'.
    Load product. Defualt format is '.dim'.
    """

    from snappy import ProductIO

    return ProductIO.readProduct(src + type)


def save_product(source, path, format='BEAM-DIMAP'):
    """
    Zapisz wynik w określonej ścieżce. Domyślny format to 'BEAM-DIMAP'.
    Save result in specified path. Defualt format is 'BEAM-DIMAP'.
    """

    from snappy import File
    from snappy import ProgressMonitor
    import jpy

    System = jpy.get_type("java.lang.System")

    PrintPM = jpy.get_type('com.bc.ceres.core.PrintWriterProgressMonitor')
    pm = PrintPM(System.out)

    WriteOp = jpy.get_type('org.esa.snap.core.gpf.common.WriteOp')
    writeOp = WriteOp(source, File(path), format)
    writeOp.writeProduct(pm)


def split_product(source, firstBurstIndex, lastBurstIndex, subswath):
    """
    Operacja split. Wycięcie obrazu SAR z określonego subswatch i indeksach serii.
    Split operation. Clip SAR image from specified subswatch and burst indexes.
    """

    from snappy import HashMap, GPF
    params = HashMap()
    params.put('firstBurstIndex', firstBurstIndex)
    params.put('lastBurstIndex', lastBurstIndex)
    params.put('subswath', subswath)

    splited = GPF.createProduct("TOPSAR-Split", params, source)
    return splited


def subset_product(source, location):
    """
    Wycięcia obszaru o zadanych współrzędnych z obrazu SAR. Współrzędne znajdują się w słowniku, przy wywoływaniu fukcji podaje się wybraną lokację ze słownika.
    Clip SAR image with specified coordinates. Coordinates are in dictionary, in calling function user specifies location drom dictionary.
    """

    from snappy import HashMap, GPF
    import jpy

    WKTReader = jpy.get_type('com.vividsolutions.jts.io.WKTReader')

    wkt_dict = {
        "wizna": 'POLYGON((22.37 53.24, 22.81 53.30, 22.81 53.54, 22.37 53.49, 22.37 53.24))',
        "sri_lanka": 'POLYGON((80.8685 7.836,81.1454 7.836,81.1454 8.0155,80.8685 8.0155,80.8685 7.836))',
        "mongolia": 'POLYGON((104.083 48.2481,103.898 47.5818,106.6222 47.2556,106.8719 48.0269,104.083 48.2481))',
        "bangladesh": 'POLYGON((90.783 24.1526,91.2714 24.0644,91.3834 24.5524,90.8784 24.6206,90.783 24.1526))'
    }
    geom = WKTReader().read(list(wkt_dict.values())[location])

    params = HashMap()
    params.put('geoRegion', geom)
    params.put('outputImageScaleInDb', False)
    params.put('copyMetadata', 'true')
    subset = GPF.createProduct("Subset", params, source)
    return subset


def subset_product_band(source, search_band):
    """
    Wybranie określonego kanału z obrazu SAR.
    Select specified band from SAR product.
    """

    from snappy import HashMap, GPF
    import jpy
    import re

    pattern = search_band
    bands = source.getBandNames()
    print("\n" + pattern)
    for band in bands:
        print(band)
        if re.match(pattern, band):
            print("\nMATCHED:" + band)
            save_band_list = band

    save_band = save_band_list
    print("\nSearched band: " + save_band + "\n")

    params = HashMap()
    params.put('sourceBands', save_band)
    subset = GPF.createProduct("Subset", params, source)
    return subset


def apply_orbit_file(source):
    """
    Pozyskiwanie pliku orbity satelitarnej, uzyskanie dokładniejszych danych o orbicie.
    Apply orbit file, get more precise orbit data.
    """

    from snappy import HashMap, GPF

    params = HashMap()
    params.put('orbitType', 'Sentinel Precise (Auto Download)')
    out = GPF.createProduct("Apply-Orbit-File", params, source)
    return out


def backgeocoding(source_master, source_slave):
    """
    Wspólna rejestracja przy użyciu daneych o orbicie i pliku DEM.
    Co-registration using orbit data and DEM file.
    """

    from snappy import HashMap, GPF

    sources = []
    sources.append(source_master)
    sources.append(source_slave)

    params = HashMap()
    params.put('demName', 'SRTM 3Sec')
    params.put('demResamplingMethod', 'BICUBIC_INTERPOLATION')
    params.put('resamplingType', 'BISINC_5_POINT_INTERPOLATION')
    out = GPF.createProduct("Back-Geocoding", params, sources)
    return out


def interferogram(source):
    """
    Oblicz interferogram z zarejestrowanych obrazów Sentinel 1.
    Compute interferogram from coregistered Sentinel 1 images.
    """

    from snappy import HashMap, GPF

    params = HashMap()
    params.put('subtractFlatEarthPhase', 'true')
    out = GPF.createProduct("Interferogram", params, source)
    return out


def get_bands(source):
    from snappy import ProductIO

    bands = source.getBandNames()
    for band in bands:
        print(band)


def terrain_correction(source, location):
    """
    Ortorektyfikacja z użyciem określonego układu współrzędych.
    Orthorectification with specified coordinate system.
    """

    from snappy import HashMap, GPF

    terrain_dict = {
        "wizna": 'EPSG:32631',
        "sri_lanka": 'EPSG:32644',
        "mongolia": 'EPSG:32647',
        "bangladesh": 'EPSG:24306'
    }

    proj = list(terrain_dict.values())[location]

    params = HashMap()
    params.put('mapProjection', proj)
    out = GPF.createProduct("Terrain-Correction", params, source)
    return out


def deburst(source):
    """
    Deburst
    """

    from snappy import HashMap, GPF

    params = HashMap()
    out = GPF.createProduct("TOPSAR-Deburst", params, source)
    return out


def dispose(sources):
    """
    Usun wyniki zapisane w pamięci.
    Remove results from memory.
    """

    for x in sources:
        x.dispose()


def calibration(source):
    """
    Kalibracja
    Calibration
    """

    from snappy import HashMap, GPF

    params = HashMap()
    out = GPF.createProduct("Calibration", params, source)
    return out


def spackle_filter(source, filter):
    """
    Filtrowanie z użyciem wybranego filtra.
    Filter with selected filter.
    """

    from snappy import HashMap, GPF

    params = HashMap()
    params.put('filter', filter)
    out = GPF.createProduct("Speckle-Filter", params, source)
    return out


def convert_datatype(source, type='uint8'):
    """
    Konwertuj typ danych.
    Convert data type.
    """

    from snappy import HashMap, GPF

    params = HashMap()
    params.put('targetDataType', type)
    out = GPF.createProduct("Convert-Datatype", params, source)
    return out


def reproject(source):
    """
    Zmień układ współrzędnych obrazu.
    Change image's coordinate system.
    """

    from snappy import HashMap, GPF

    params = HashMap()
    # params.put('mapProjection', 'EPSG:32631 - WGS 84 / UTM zone 31N')
    params.put('crs', 'EPSG:32631')
    out = GPF.createProduct("Reproject", params, source)
    return out


def resampling(source):
    from snappy import HashMap, GPF

    params = HashMap()
    out = GPF.createProduct("Resampling", params, source)
    return out


##############################################################

def read_product_cv2(source, ext='.tif'):
    """
    Wczytaj obraz używając openCV w skali szarości
    Load image using openCV in grayscale
    """

    import cv2

    source += ext
    img = cv2.imread(source, 0)
    print(source)
    return img


def save_product_cv2(source, path):
    """
    Zapisz obraz używając openCV.
    Save image using openCV.
    """

    import cv2

    path += ".tif"

    cv2.imwrite(path, source)


def save_product_grass(source, path):
    """
    Zapisz obraz używając GRASS GIS.
    Save product using GRASS GIS.
    """

    import grass_session as session
    import grass.script as gscript

    path += '.tif'
    gscript.run_command('r.in.gdal', input=source, output='product', flags='o', overwrite=True)
    gscript.run_command('r.out.gdal', input='product', output='Grass.tif', overwrite=True)


def thresholding_otsu(source):
    """
    Progowanie metodą Otsu.
    Thresholding with Otsu method.
    """

    import cv2

    blur = cv2.GaussianBlur(source, (5, 5), 0)
    ret3, th3 = cv2.threshold(blur, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    return th3


def projection(dataset_source, dataset_out):
    """
    Pobierz układ współrzędych jednego obrazu i zastosuj go w drugim obrazie.
    Get coordinate system from one image and apply it to onether one.
    """
    import gdal

    dataset1 = gdal.Open(dataset_source)
    projection = dataset1.GetProjection()
    geotransform = dataset1.GetGeoTransform()

    dataset2 = gdal.Open(dataset_out, gdal.GA_Update)
    dataset2.SetGeoTransform(geotransform)
    dataset2.SetProjection(projection)


def hand_calculator(source, path, project, mapset):
    """
    Oblicz HAND
    Calculate HAND
    """

    from pysheds.grid import Grid
    import grass_session as session
    import grass.script as gscript
    import gdal
    from gdalconst import GA_Update
    import shutil

    path_grass = 'C:\\Users\\micha\\Documents\\grassdata\\'
    # read
    try:
        gscript.run_command('g.proj', georef=source, flags='ct', location=mapset)
        print("\nLocation " + mapset + " from " + source + " created.\n")
    except:
        print(f"Location already exsists. Remove location {mapset}.")
        shutil.rmtree(path_grass + mapset)
        print(f"Location {mapset} removed.")
        gscript.run_command('g.proj', georef=source, flags='ct', location=mapset)
        print("\nLocation " + mapset + " from " + source + " created.\n")

    gscript.run_command('g.mapset', mapset='PERMANENT', location=mapset)
    gscript.run_command('r.in.gdal', input=source, output='dem', flags='o', overwrite=True)
    gscript.run_command('r.out.gdal', input='dem', output='dem.tif', overwrite=True)

    thr = 70000
    thr_loop = 10000
    lim = 160000
    loop = False
    if (loop):
        while (thr_loop < lim):
            # stream
            thr = thr_loop
            out_stream = path + '\\' + project + '_stream' + str(thr) + '.tif'
            gscript.run_command('r.watershed', elevation='dem', threshold=thr, drainage='dem_drain', basin='dem_basin',
                                stream='dem_stream', overwrite=True)
            gscript.run_command('r.null', map='dem_stream', null=0)
            gscript.run_command('r.out.gdal', input='dem_stream', output=out_stream, overwrite=True)
            projection(source, out_stream)
            # hand
            grid = Grid.from_raster(source, data_name='dem')
            # Fill depressions in DEM
            grid.fill_depressions('dem', out_name='flooded_dem')
            # Resolve flats in DEM
            grid.resolve_flats('flooded_dem', out_name='inflated_dem')
            # Specify directional mapping
            dirmap = (3.0, 2.0, 1.0, 8.0, 7.0, 6.0, 5.0, 4.0)
            # Compute flow directions
            grid.flowdir(data='inflated_dem', out_name='dir', dirmap=dirmap)
            grid.read_raster(out_stream, data_name='drainage')
            grid.compute_hand(fdir='dir', dem='dem', drainage_mask='drainage', out_name='hand', dirmap=dirmap,
                              nodata_in_fdir=0)
            out_hand = f"{path}\\{project}result_{str(thr)}.tif"
            grid.to_raster('hand', out_hand)

            thr_loop += 10000
            ds = gdal.Open(out_hand, GA_Update)
            dsBand = ds.GetRasterBand(1)
            result = gdal.FillNodata(targetBand=dsBand, maskBand=None,
                                     maxSearchDist=5, smoothingIterations=0)
            del ds
    # stream
    out_stream = path + '\\' + project + '_stream' + str(thr) + '.tif'
    gscript.run_command('r.watershed', elevation='dem', threshold=thr, drainage='dem_drain', basin='dem_basin',
                        stream='dem_stream', overwrite=True)
    gscript.run_command('r.null', map='dem_stream', null=0)
    gscript.run_command('r.out.gdal', input='dem_stream', output=out_stream, overwrite=True)
    projection(source, out_stream)
    # hand
    grid = Grid.from_raster(source, data_name='dem')
    # Fill depressions in DEM
    grid.fill_depressions('dem', out_name='flooded_dem')
    # Resolve flats in DEM
    grid.resolve_flats('flooded_dem', out_name='inflated_dem')
    # Specify directional mapping
    dirmap = (3.0, 2.0, 1.0, 8.0, 7.0, 6.0, 5.0, 4.0)
    # Compute flow directions
    grid.flowdir(data='inflated_dem', out_name='dir', dirmap=dirmap)
    grid.read_raster(out_stream, data_name='drainage')
    grid.compute_hand(fdir='dir', dem='dem', drainage_mask='drainage', out_name='hand', dirmap=dirmap)
    out_hand = f"{path}\\{project}_result.tif"
    grid.to_raster('hand', out_hand)
    ds = gdal.Open(out_hand, GA_Update)
    dsBand = ds.GetRasterBand(1)
    result = gdal.FillNodata(targetBand=dsBand, maskBand=None,
                             maxSearchDist=5, smoothingIterations=0)
    del ds


def map_calculator(files, hand_file, outpath, mapset):
    """
    Maskowanie obrazu.
    Mask images.
    """

    import gdal
    import grass_session as session
    import grass.script as gscript
    import os
    import shutil

    if not os.path.exists(outpath):
        os.mkdir(outpath)

    path_grass = 'C:\\Users\\micha\\Documents\\grassdata\\'
    try:
        gscript.run_command('g.proj', georef=hand_file, flags='ct', location=mapset)
        print("\nLocation " + mapset + " created.\n")
    except:
        print(f"Location already exsists. Remove location {mapset}.")
        shutil.rmtree(path_grass + mapset)
        print(f"Location {mapset} removed.")
        gscript.run_command('g.proj', georef=hand_file, flags='ct', location=mapset)
        print("\nLocation " + mapset + " created.\n")

    gscript.run_command('g.mapset', mapset='PERMANENT', location=mapset)
    input_dict = {}
    expressions = {
        # dla sri lanki hand <= 20
        # dla mongolii grd_vh <= 0.005, grd_vh <=0.015
        'otsu': 'result = if("source" == 0 && "mask_hand" <= 5)',
        'pred': 'result = if("source" == 1 && "mask_hand" <= 5)',
        'slc_interferogram': 'result = if(("source" <= 0.25 && "source" > 0) && "mask_hand" <= 5)',
        'grd_vh': 'result = if(("source" <= 0.015 && "source" > 0) && "mask_hand" <= 5)',
        'grd_vv': 'result = if(("source" <= 0.06 && "source" > 0) && "mask_hand" <= 5)'
        # 'grd_vh':'result = if(("source" <= 0.005 && "source" != 0) && "mask_hand" <= 5)',
        # 'grd_vv':'result = if(("source" <= 0.015 && "source" != 0) && "mask_hand" <= 5)'
    }
    for f in files:
        f_low = f.lower()
        if 'otsu' in f_low:
            input_dict[f] = expressions['otsu']
        if 'pred' in f_low:
            input_dict[f] = expressions['pred']
        if 'interferogram' in f_low and 'otsu' not in f_low:
            input_dict[f] = expressions['slc_interferogram']
        if (('grd' in f_low and 'otsu' not in f_low) and ('vh' in f_low)):
            input_dict[f] = expressions['grd_vh']
        if (('grd' in f_low and 'otsu' not in f_low) and ('vv' in f_low)):
            input_dict[f] = expressions['grd_vv']
    for f in files:
        sep = "\\"
        out = f"{outpath}{f.split(sep)[-1][:-4]}_Result.tif"
        print('\nRASTER: ' + f)
        print('EXPRESSION: ' + input_dict[f])
        gscript.run_command('r.in.gdal', input=f, output='source', flags='o', overwrite=True)
        gscript.run_command('g.region', raster='source')
        gscript.run_command('r.in.gdal', input=hand_file, output='mask_hand', flags='o', overwrite=True)
        gscript.run_command('r.mapcalc', expression=input_dict[f], overwrite=True)
        gscript.run_command('r.null', map='result', null=0)
        gscript.run_command('r.out.gdal', input='result', output=out, overwrite=True)
        projection(f, out)
