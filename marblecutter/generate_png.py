import gdal, gdalconst
import numpy as np
import sys
import math
import os

def compute_cutoff(indata, percentage, enviLike=False):
    data = None
    if enviLike:
        # data = data[ int(data.shape[0]/2-500):int(data.shape[0]/2+500), int(data.shape[1]/2-500):int(data.shape[1]/2+1000) ]
        data = indata[::100,::100]
    else:
        data = indata
    hist = list(range(256))
    minMax = ( data[data>0].min(), data[data>0].max())
    stepSize = (minMax[1] - minMax[0])/256.0
    for i in range(255):
        binMin = minMax[0] + stepSize*i
        binMax = minMax[0] + stepSize*(i+1)
        hist[i] = np.sum( data[data>=binMin] < binMax )
    hist[255] = np.sum( data >= (minMax[0] + stepSize*255 ) )
    minCutoff = minMax[0]
    maxCutoff = minMax[1]
    totalCount = np.float(np.sum( hist ))
    cumHist = list(range(256))
    cumHist[0] = hist[0]
    minCalculated = False
    maxCalculated = False
    for i in range(1,256):
        cumHist[i] = cumHist[i-1] + hist[i]
        if np.float(cumHist[i]/totalCount)*100.0 >= np.float(percentage) and minCalculated==False:
            print("HEY: ", np.float(cumHist[i]/totalCount)*100.0)
            minCutoff = minMax[0] + stepSize*i
            minCalculated = True
        if np.float(cumHist[i]/totalCount)*100.0 >= np.float(100-percentage) and maxCalculated==False:
            maxCutoff = minMax[0] + stepSize*i
            maxCalculated = True
    print( minCutoff, maxCutoff)    
    return (minCutoff, maxCutoff)

def generate_8bit_3857(bandFilePaths,  productCode, mergeFileDir):
    ds = gdal.Open(bandFilePaths[0])
    dataType = ds.GetRasterBand(1).DataType
    numBands = len(bandFilePaths)
    driver = gdal.GetDriverByName("GTiff")
    mergeTifFilename = os.path.join(mergeFileDir, "{0}_OSR_merge.tif".format(productCode))
    merge8bitTifFilename = os.path.join(mergeFileDir, "{0}_OSR_merge_8bit.tif".format(productCode)) 
    merge8bit3857TifFilename = os.path.join(mergeFileDir, "{0}_OSR_merge_8bit_EPSG3857.tif".format(productCode))
    tmpDs = driver.Create(mergeTifFilename, ds.RasterXSize, ds.RasterYSize, numBands, dataType)
    bandWiseMinMax = list(range(numBands))
    for i in range(numBands):
        inDs = gdal.Open(bandFilePaths[i])
        inData = inDs.GetRasterBand(1).ReadAsArray()
        bandWiseMinMax[i] = compute_cutoff(inData, 2, True)
        tmpDs.GetRasterBand(i+1).WriteArray(inData)
    tmpDs.SetProjection(ds.GetProjection())
    tmpDs.SetGeoTransform(ds.GetGeoTransform())
    tmpDs = None
    cmd = "gdal_translate -ot Byte "
    for i in range(numBands):
        cmd += "-scale_{0} {1} {2} 1 254 ".format(i+1, int(bandWiseMinMax[i][0]) , int(bandWiseMinMax[i][1]))
    cmd += "-of GTiff {0} {1}".format(mergeTifFilename, merge8bitTifFilename)
    os.system(cmd) 
    if os.path.exists(merge8bit3857TifFilename):
        os.unlink(merge8bit3857TifFilename)
    cmd = "gdalwarp -of GTiff -t_srs \"EPSG:3857\" -srcnodata \""
    for i in range(numBands):
        cmd += "0 "
    cmd += "\" -dstnodata \""   
    for i in range(numBands):
        cmd += "0 "
    cmd += "\" {0} {1}".format(merge8bitTifFilename, merge8bit3857TifFilename)
    os.system(cmd)
    return merge8bitTifFilename

def generate_png(tifFilePath, pngFilePath):
    cmd = "gdal_translate -of PNG -outsize 512 512 -b 3 -b 2 -b 1 {0} {1}".format(tifFilePath, pngFilePath)
    os.system(cmd)
    
def generate_png_back(band2FilePath, band3FilePath, band4FilePath, pngFilePath):
    driver    = gdal.GetDriverByName("GTiff")
    band2Ds   = gdal.Open(band2FilePath)
    band3Ds   = gdal.Open(band3FilePath)
    band4Ds   = gdal.Open(band4FilePath)
    outsize   = 512
    inXSize    = band2Ds.RasterXSize
    inYSize    = band2Ds.RasterYSize
    stepXSize  = math.ceil(inXSize/outsize)
    stepYSize  = math.ceil(inYSize/outsize)
    tifFilePath = pngFilePath[:-3]+"tif"
    alphaTif    = pngFilePath[:-4]+"_alpha.tif"
    outDs     = driver.Create(tifFilePath, outsize, outsize, 3, gdalconst.GDT_Byte)
    inBand    = [ band4Ds.GetRasterBand(1) ,
        band3Ds.GetRasterBand(1),
        band2Ds.GetRasterBand(1) ]
    for i in range(3):
        outBand = outDs.GetRasterBand(i+1)
        outData = np.zeros((512, 512))
        inData = inBand[i].ReadAsArray()[::stepYSize, ::stepXSize]
        outData[:inData.shape[0], :inData.shape[1]] = np.uint8 ( (inData - inData.min())/(inData.max() - inData.min()) * 345)
        outBand.WriteArray(outData)
    outDs.SetProjection(band2Ds.GetProjection())
    outDs.SetGeoTransform(band2Ds.GetGeoTransform())
    outDs = None
    command = "gdal_translate -of PNG " + tifFilePath + " " + pngFilePath
    
    os.system(command)
    os.unlink(tifFilePath)

    return True

if __name__ == '__main__':
    band2FilePath = sys.argv[1]
    band3FilePath = sys.argv[2]
    band4FilePath = sys.argv[3]
    band5FilePath = sys.argv[4]
    pngFilePath = sys.argv[5]
    mergeFileDir = sys.argv[6]
    productCode = sys.argv[7]
    # generate_png(band2FilePath, band3FilePath, band4FilePath, pngFilePath)
    tifFilePath = generate_8bit_3857([band2FilePath,band3FilePath,band4FilePath,band5FilePath], productCode, mergeFileDir)
    generate_png(tifFilePath, pngFilePath)
