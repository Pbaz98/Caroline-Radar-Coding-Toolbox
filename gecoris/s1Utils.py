#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""s1Utils

Module containing Sentinel-1 utility functions

Copyright (C) 2020 by R.Czikhardt

Email: czikhardt.richard@gmail.com
Last edit: 28.7.2020

This file is part of GECORIS - Geodetic Corner Reflector (In)SAR Toolbox.

    GECORIS is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GECORIS is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with GECORIS. If not, see <https://www.gnu.org/licenses/>.
"""

import numpy as np
# from snappy import jpy
# from snappy import ProductIO
from datetime import datetime
# gecoris packages:
from gecoris import geoUtils, radarUtils
# constants:
speedOfLight = 299792458.0


def readMetadata(file,MS = 'master'):
    """Return a metadata dictionary
    
    Evaluate orbit state vector at given azimuth time using fitted Chebyshev 
    polynomials by function 'orbitFit'.
    
    input: 'orbitFit' dict (as generated by 'orbitFit' function)
            MS - 'master'/'slave'
    """    
    
    # get necesarry Java-Python models:
    AbstractMetadata = jpy.get_type('org.esa.snap.engine_utilities.datamodel.AbstractMetadata')

    # read coregistered product:
    product = ProductIO.readProduct(file)
    
    # read metadata root:
    if MS == 'master':
        AbsRoot = AbstractMetadata.getAbstractedMetadata(product)
    elif MS == 'slave':
        AbsRoot = AbstractMetadata.getSlaveMetadata(product.getMetadataRoot()).getElementAt(0)
    else:
        raise('Must specify either "master" or "slave".')
    
    # get required metadata:
    
    #acqDate  = AbsRoot.getAttributeString(AbstractMetadata.first_line_time)
    #acqDate = datetime.strptime(acqDate[0:20], '%d-%b-%Y %H:%M:%S') - problem with locale
    acqDate = AbsRoot.getAttributeString('PRODUCT')[17:32]
    acqDate = datetime.strptime(acqDate,'%Y%m%dT%H%M%S')
    
    geometry = AbsRoot.getAttributeString('PASS')
    
    firstLineTime = AbstractMetadata.parseUTC(AbsRoot.getAttributeString(AbstractMetadata.first_line_time)).getMJD()
    # convert to [seconds] of day:
    firstLineTime = np.mod(firstLineTime,1)*24*3600
    
    range0time = AbsRoot.getAttributeDouble(AbstractMetadata.slant_range_to_first_pixel)/speedOfLight
    
    PRF = 1/AbsRoot.getAttributeDouble(AbstractMetadata.line_time_interval)
    
    RSR = AbsRoot.getAttributeDouble(AbstractMetadata.range_sampling_rate)*2*1e6
    
    wavelength = speedOfLight/(AbsRoot.getAttributeDouble(AbstractMetadata.radar_frequency)*1e6)
    
    rangeSampling = AbsRoot.getAttributeDouble(AbstractMetadata.range_spacing)
    
    azimuthSampling = AbsRoot.getAttributeDouble(AbstractMetadata.azimuth_spacing)
    
    nAzimuth = AbsRoot.getAttributeInt(AbstractMetadata.num_output_lines)
    nRange = AbsRoot.getAttributeInt(AbstractMetadata.num_samples_per_line)
    
    # centre coordinates:    
    centerLat = (AbsRoot.getAttributeDouble(AbstractMetadata.first_near_lat) +
                 AbsRoot.getAttributeDouble(AbstractMetadata.first_far_lat) +
                 AbsRoot.getAttributeDouble(AbstractMetadata.last_near_lat) +
                 AbsRoot.getAttributeDouble(AbstractMetadata.last_far_lat))/4
    centerLon = (AbsRoot.getAttributeDouble(AbstractMetadata.first_near_long) +
                 AbsRoot.getAttributeDouble(AbstractMetadata.first_far_long) +
                 AbsRoot.getAttributeDouble(AbstractMetadata.last_near_long) +
                 AbsRoot.getAttributeDouble(AbstractMetadata.last_far_long))/4
    centerH = AbsRoot.getAttributeDouble(AbstractMetadata.avg_scene_height)
    centerAzimuth = np.round(nAzimuth/2)
    
    swath = int(AbsRoot.getAttributeString(AbstractMetadata.SWATH)[-1])
       
    # get Orbit:
    orbit = AbstractMetadata.getOrbitStateVectors(AbsRoot)
    # fit orbit data:
    orbFit = geoUtils.orbitFit(orbit,verbose = 0)
    
    # get beta0 calibration factor:
    beta0 = getBeta0Cal(file)
    
    # resolutions [from s1 annual performance reports]
    azResolutions = np.array([21.76,21.89,21.71])
    rResolutions = np.array([2.63,3.09,3.51])
    
    s1utils = jpy.get_type('org.esa.s1tbx.commons.Sentinel1Utils')
    su = s1utils(product)
    subSwath = su.getSubSwath()
    # chcek if not multi-swath product:
    if subSwath:
        nBursts = subSwath[0].numOfBursts
        # get burst metadata if raw product:
        if nBursts > 1:
    
            firstAzTime = np.array(subSwath[0].burstFirstLineTime)
            firstAzTime = np.mod(firstAzTime/86400,1)*24*3600
            firstRgTime = subSwath[0].slrTimeToFirstPixel
            nAzBurst = subSwath[0].linesPerBurst
            nRgBurst = subSwath[0].samplesPerBurst
            nAzTotal = subSwath[0].numOfLines
            nRgTotal = subSwath[0].numOfSamples
            # burst metadata dict:
            burstInfo = {
                    "nBursts": nBursts,
                    "firstAzTime": firstAzTime,
                    "firstRgTime": firstRgTime,
                    "nAzBurst": nAzBurst,
                    "nRgBurst": nRgBurst,
                    "nAzTotal": nAzTotal,
                    "nRgTotal": nRgTotal}
        else:
             burstInfo = None
        # for S1 get TOPS metadata:
        # 1.) steering rate
        steeringRate = subSwath[0].azimuthSteeringRate*np.pi/180
        
        azRes = azResolutions[swath-1]
        rRes = rResolutions[swath-1]
    else:
        burstInfo = None
        steeringRate = None
        nBursts = 0
        swath = 'ALL'
        azRes = azResolutions[1]
        rRes = rResolutions[1]
        
    # for S1 get TOPS metadata:
    # 2.) get azimuth fm rate polynomial for all bursts:
    azFmRateArray = getAzFmRate(product)
    # 3.) get Doppler centroid:
    dcPolyArray = getDcPoly(product)
    
    # get downlink info:
    (pri,rank,chirpRate) = getDownlinkInfo(product)
    
    # Create metadata dictionary:
    metadata = {
      "orbit": geometry,
      "acqDate": acqDate,
      "azimuth0time": firstLineTime,
      "range0time": range0time,
      "PRF": PRF,
      "RSR": RSR,
      "wavelength": wavelength,
      "orbitFit": orbFit,
      "rangeSpacing": rangeSampling,
      "azimuthSpacing": azimuthSampling,
      "centerLon": centerLon,
      "centerLat": centerLat,
      "centerH": centerH,
      "nAzimuth": nAzimuth,
      "nRange": nRange,
      "swath": swath,
      "centerAzimuth": centerAzimuth,
      "beta0": beta0,
      "azimuthResolution": azRes,
      "rangeResolution": rRes,
      "nBursts": nBursts,
      "burstInfo": burstInfo,
      "steeringRate": steeringRate,
      "azFmRateArray": azFmRateArray,
      "dcPolyArray": dcPolyArray,
      "PRI": pri,
      "rank": rank,
      "chirpRate": chirpRate
    }
    product.dispose()
    return metadata


def getBurstIdx(metadata,plh):
    """Return burst index of stack for given elipsoidal coordinates.
    
    input:
    """   
    if metadata['burstInfo']:
        xyz = geoUtils.plh2xyz(plh)
        tAz,_,_ = radarUtils.xyz2t(xyz,metadata)
        burstAz = metadata['burstInfo']['firstAzTime']
        burstIdx = np.nonzero(tAz > burstAz)[0][-1]
        if burstIdx >= 0:
            return burstIdx
        else:
            print('Point in none of the bursts')
            return
    else:
        print('No burst-wise stack!')
        return
    

def getBeta0Cal(file):
    """Return a beta0 calibration factor
    
    Read single beta0 calibration factor from SNAP SLC/coregistered product
    at 'file'. FUnction assumes homogenous beta0 over whole Sentinel-1 swath.
    
    input: 'file' full path to SNAP .dim file
    """   
    
    product = ProductIO.readProduct(file)    
    AbstractMetadata = jpy.get_type('org.esa.snap.engine_utilities.datamodel.AbstractMetadata')
    
    origMeta = AbstractMetadata.getOriginalProductMetadata(product)
    srcCalibration = origMeta.getElement("calibration")
    bandCalibration = srcCalibration.getElementAt(0).getElement("calibration");
    calibrationVectorListElem = bandCalibration.getElement("calibrationVectorList");
    calList = calibrationVectorListElem.getElements();
    # sufficient to use only first entry as beta0 doesnt change per epoch
    beta0 = list(calList)[0].getElement('betaNought').getAttributeString('betaNought')
    beta0 = float(beta0.split(' ')[0])
    
    product.dispose()
    return beta0


def getAzFmRate(product):
    #% get azimuth fm rate

    AbstractMetadata = jpy.get_type('org.esa.snap.engine_utilities.datamodel.AbstractMetadata')
    s1utils = jpy.get_type('org.esa.s1tbx.commons.Sentinel1Utils')

    origProdRoot = AbstractMetadata.getOriginalProductMetadata(product);
    annotation = origProdRoot.getElement("annotation")
    elem = annotation.getElements()[0]
    prod = elem.getElement("product");
    generalAnnotation = prod.getElement("generalAnnotation");
    azimuthFmRateList = generalAnnotation.getElement("azimuthFmRateList");
    #count = int(azimuthFmRateList.getAttributeString("count"))
    
    azFmRateListElem = azimuthFmRateList.getElements()
    azFmRateArray = [] # loop through all bursts:
    for elem in azFmRateListElem:
        azFmPoly = elem.getElement("azimuthFmRatePolynomial")
        coeffs = azFmPoly.getAttributeString("azimuthFmRatePolynomial").split(' ')
        azFmRate = {
            't' : s1utils.getTime(elem,"azimuthTime").getMJD()*86400,
            't0' : float(elem.getAttributeString("t0")),
            'c0' : float(coeffs[0]),
            'c1' : float(coeffs[1]),
            'c2' : float(coeffs[2])}
        azFmRateArray.append(azFmRate)
 
    return azFmRateArray


def getDcPoly(product):
    #% get Doppler-centroid polynomials

    AbstractMetadata = jpy.get_type('org.esa.snap.engine_utilities.datamodel.AbstractMetadata')
    s1utils = jpy.get_type('org.esa.s1tbx.commons.Sentinel1Utils')

    origProdRoot = AbstractMetadata.getOriginalProductMetadata(product)
    annotation = origProdRoot.getElement("annotation")
    elem = annotation.getElements()[0]
    prod = elem.getElement("product");
    imageAnnotation = prod.getElement("imageAnnotation");
    procInfo = imageAnnotation.getElement("processingInformation")
    dcMethod = procInfo.getAttributeString("dcMethod")
    dopplerCentroid = prod.getElement("dopplerCentroid")
    dcEstimateList = dopplerCentroid.getElement("dcEstimateList")
    #count = int(dcEstimateList.getAttributeString("count"))
    
    dcEstimateListElem = dcEstimateList.getElements()
    dcPolyArray = [] # loop through all bursts:
    for elem in dcEstimateListElem:
        t = s1utils.getTime(elem,"azimuthTime").getMJD()*86400
        t0 = elem.getAttributeDouble("t0")
        if dcMethod == 'Data Analysis':
            dcPoly = elem.getElement("dataDcPolynomial")
            coeffs = dcPoly.getAttributeString("dataDcPolynomial").split(' ')
        else:
            dcPoly = elem.getElement("geometryDcPolynomial")
            coeffs = dcPoly.getAttributeString("geometryDcPolynomial").split(' ')
        dcPoly = {
            't' : t,
            't0' : t0,
            'c0' : float(coeffs[0]),
            'c1' : float(coeffs[1]),
            'c2' : float(coeffs[2])}
        dcPolyArray.append(dcPoly)
       
    return dcPolyArray


def getDownlinkInfo(product):
    #% get downlink info

    AbstractMetadata = jpy.get_type('org.esa.snap.engine_utilities.datamodel.AbstractMetadata')
    origProdRoot = AbstractMetadata.getOriginalProductMetadata(product)
    annotation = origProdRoot.getElement("annotation")
    elem = annotation.getElements()[0]
    prod = elem.getElement("product");
    generalAnnotation = prod.getElement("generalAnnotation");
    
    downlinkList = generalAnnotation.getElement("downlinkInformationList")
    downlinkInfo = downlinkList.getElement('downlinkInformation')
    downlinkValues = downlinkInfo.getElement('downlinkValues')
    pri = float(downlinkValues.getAttributeString("pri"))
    rank = float(downlinkValues.getAttributeString("rank"))
    chirpRate = float(downlinkValues.getAttributeString("txPulseRampRate"))
      
    return (pri,rank,chirpRate)


def getDerampDemodPhase(metadata,boundingBox,ovsFactor = 1):
    
    # parse input:    
    minAz = boundingBox[0][0]
    maxAz = boundingBox[0][1]
    minR = boundingBox[1][0]
    maxR = boundingBox[1][1]
    
    # prepare time vectors:
    rangeVec = np.arange(minR,maxR+1)
    azimuthVec = np.arange(minAz,maxAz+1)    
    (tAzimuth,tRange) = radarUtils.radar2time(azimuthVec,rangeVec,metadata)
    
    # OVS:
    if ovsFactor > 1:
        nodes = np.arange(0,len(tAzimuth))*ovsFactor#+ovsFactor/2
        ovsNodes = np.arange(0,(len(tAzimuth)*ovsFactor))
        tAzimuth = np.polyval(np.polyfit(nodes, tAzimuth, 1),ovsNodes)
        nodes = np.arange(0,len(tRange))*ovsFactor#+ovsFactor/2
        ovsNodes = np.arange(0,(len(tRange)*ovsFactor))
        tRange = np.polyval(np.polyfit(nodes, tRange, 1),ovsNodes)
    
    # correct burst timing:
    if metadata['nBursts'] > 1:
        burstAz = metadata['burstInfo']['firstAzTime']
        idx = np.nonzero(np.mean(tAzimuth) > burstAz)[0][-1]
        tmp1 = metadata['burstInfo']['nAzBurst']/2/metadata['PRF'] # burst +-time span (Eta)
        midAzimuthTime = burstAz[idx] + tmp1 # burst +-time span (Eta)
    else:
        tmp1 = metadata['nAzimuth']/2/metadata['PRF'] # burst +-time span (Eta)
        midAzimuthTime = metadata['azimuth0time'] + tmp1 # burst +-time span (Eta)
        
    # velocity:
    dummy,vel = geoUtils.orbitVal(metadata['orbitFit'],midAzimuthTime)
    vel = np.linalg.norm(vel) # absolute value
    # angular steering rate:
    krot = 2*vel*metadata['steeringRate']/metadata['wavelength']
    
    # get azimuth fm rate for all bursts:
    azFmRateArray = metadata['azFmRateArray']
    # burst central times:
    tBurst = [np.mod(tmp['t']/3600/24,1)*3600*24 for tmp in azFmRateArray]
    # find appropriate burst index: 
    burstIdx = np.argmin(np.abs(tBurst[0]-midAzimuthTime))
    #print('WARNING: addad [0] to tBurst !! to remove or adjust because before tBurst was a list. why??')
    
    # compute range dependent DopplerRate:
    dt = tRange*2 - azFmRateArray[burstIdx]['t0']
    rangeDependDopplerRate = (azFmRateArray[burstIdx]['c0'] + 
                              azFmRateArray[burstIdx]['c1']*dt + 
                              azFmRateArray[burstIdx]['c2']*dt*dt)
    # DopperRate [Hz/s]:
    dopplerRate = rangeDependDopplerRate*krot/(rangeDependDopplerRate - krot)
    
    #% get Doppler centroid:
    #dcPolyArray = getDcPoly(product)
    dcPolyArray = metadata['dcPolyArray']
    dt = tRange*2 - dcPolyArray[burstIdx]['t0']
    dopplerCentroid = (dcPolyArray[burstIdx]['c0'] + 
                       dcPolyArray[burstIdx]['c1']*dt + 
                       dcPolyArray[burstIdx]['c2']*dt*dt)
    
    # Doppler centroid and FM-rate at first range sample of burst:
    ## reference time computation:
    midRangeTime = metadata['range0time']*2+metadata['nRange']/2/metadata['RSR']
    
    dt0 = midRangeTime - dcPolyArray[burstIdx]['t0']
    dc0 = (dcPolyArray[burstIdx]['c0'] + dcPolyArray[burstIdx]['c1']*dt0 + 
       dcPolyArray[burstIdx]['c2']*dt0*dt0)
    dt0 = midRangeTime - azFmRateArray[burstIdx]['t0']
    fm0 = (azFmRateArray[burstIdx]['c0'] + azFmRateArray[burstIdx]['c1']*dt0 + 
       azFmRateArray[burstIdx]['c2']*dt0*dt0)
    # finally reference time:
    refTime = dopplerCentroid/rangeDependDopplerRate- dc0/fm0
    
    # finally compute deramp phase
    ## prepare grids for vectorized computation:    
    azGrid = np.tile(tAzimuth,(tRange.size,1)).T - midAzimuthTime
    refTimeGrid = np.tile(refTime,(tAzimuth.size,1))
    dopplerRateGrid =  np.tile(dopplerRate,(tAzimuth.size,1))
    dopplerCentroidGrid = np.tile(dopplerCentroid,(tAzimuth.size,1))
    
    derampPhase = -np.pi * dopplerRateGrid * np.power(azGrid - refTimeGrid,2)
    demodPhase = -2*np.pi * dopplerCentroidGrid*(azGrid - refTimeGrid)
    
    return derampPhase+demodPhase


def computeDerampDemodPhase(metadata,boundingBox,ovsFactor = 1):
    
    # parse input:
    minAz = boundingBox[0][0]
    maxAz = boundingBox[0][1]
    minR = boundingBox[1][0]
    maxR = boundingBox[1][1]
    
    # prepare time vectors:
    rangeVec = np.arange(minR,maxR+1)
    azimuthVec = np.arange(minAz,maxAz+1)    
    (tAzimuth,tRange) = radarUtils.radar2time(azimuthVec,rangeVec,metadata)
    
    # OVS:
    if ovsFactor > 1:
        nodes = np.arange(0,len(tAzimuth))*ovsFactor#+ovsFactor/2
        ovsNodes = np.arange(0,(len(tAzimuth)*ovsFactor))
        tAzimuth = np.polyval(np.polyfit(nodes, tAzimuth, 1),ovsNodes)
        nodes = np.arange(0,len(tRange))*ovsFactor#+ovsFactor/2
        ovsNodes = np.arange(0,(len(tRange)*ovsFactor))
        tRange = np.polyval(np.polyfit(nodes, tRange, 1),ovsNodes) 
 
    # mid-burst azimuth time
    if metadata['nBursts'] > 1:
        burstAz = metadata['burstInfo']['firstAzTime']
        idx = np.nonzero(np.mean(tAzimuth) > burstAz)[0][-1]
        tmp1 = metadata['burstInfo']['nAzBurst']/2/metadata['PRF'] # burst +-time span (Eta)
        midAzimuthTime = burstAz[idx] + tmp1 # burst +-time span (Eta)
    else:
        tmp1 = metadata['nAzimuth']/2/metadata['PRF'] # burst +-time span (Eta)
        midAzimuthTime = metadata['azimuth0time'] + tmp1 # burst +-time span (Eta)
    
    # velocity:
    dummy,vel = geoUtils.orbitVal(metadata['orbitFit'],midAzimuthTime)
    vel = np.linalg.norm(vel) # absolute value
    # angular steering rate:
    krot = 2*vel*metadata['steeringRate']/metadata['wavelength']
    
    # azimuth fm rate for all bursts:
    azFmRateArray = metadata['azFmRateArray']
    # burst central times:
    tBurst = [np.mod(tmp['t']/3600/24,1)*3600*24 for tmp in azFmRateArray]
    # find appropriate burst index: 
    burstIdx = np.argmin(np.abs(tBurst-midAzimuthTime))
    
    # compute range dependent DopplerRate:
    dt = tRange*2 - azFmRateArray[burstIdx]['t0']
    rangeDependDopplerRate = (azFmRateArray[burstIdx]['c0'] + 
                              azFmRateArray[burstIdx]['c1']*dt + 
                              azFmRateArray[burstIdx]['c2']*dt*dt)
    # DopperRate [Hz/s]:
    dopplerRate = rangeDependDopplerRate*krot/(rangeDependDopplerRate - krot)
    
    #% get Doppler centroid:
    dcPolyArray = metadata['dcPolyArray']
    dt = tRange*2 - dcPolyArray[burstIdx]['t0']
    dopplerCentroid = (dcPolyArray[burstIdx]['c0'] + 
                       dcPolyArray[burstIdx]['c1']*dt + 
                       dcPolyArray[burstIdx]['c2']*dt*dt)
    
    # Doppler centroid and FM-rate at first range sample of burst:
    midRangeTime = metadata['range0time']*2+metadata['nRange']/2/metadata['RSR']    
    dt0 = midRangeTime - dcPolyArray[burstIdx]['t0']
    dc0 = (dcPolyArray[burstIdx]['c0'] + dcPolyArray[burstIdx]['c1']*dt0 + 
       dcPolyArray[burstIdx]['c2']*dt0*dt0)
    dt0 = midRangeTime - azFmRateArray[burstIdx]['t0']
    fm0 = (azFmRateArray[burstIdx]['c0'] + azFmRateArray[burstIdx]['c1']*dt0 + 
       azFmRateArray[burstIdx]['c2']*dt0*dt0)
    # finally reference time:
    refTime = dopplerCentroid/rangeDependDopplerRate- dc0/fm0
    
    # finally compute deramp phase
    ## prepare grids for vectorized computation:    
    azGrid = np.tile(tAzimuth,(tRange.size,1)).T - midAzimuthTime
    refTimeGrid = np.tile(refTime,(tAzimuth.size,1))
    dopplerRateGrid =  np.tile(dopplerRate,(tAzimuth.size,1))
    dopplerCentroidGrid = np.tile(dopplerCentroid,(tAzimuth.size,1))
    
    derampPhase = -np.pi * dopplerRateGrid * np.power(azGrid - refTimeGrid,2)
    demodPhase = -2*np.pi * dopplerCentroidGrid*(azGrid - refTimeGrid)
    
    return derampPhase+demodPhase


def reramp(SLC,dPhase):    
    bandI = np.real(SLC)
    bandQ = np.imag(SLC)
    I = bandI*np.cos(dPhase) + bandQ*np.sin(dPhase)
    Q = -bandI*np.sin(dPhase) + bandQ*np.cos(dPhase)
    SLCreramp = I + 1j*Q    
    return SLCreramp


def deramp(SLC,dPhase):
    bandI = np.real(SLC)
    bandQ = np.imag(SLC)    
    I = bandI*np.cos(dPhase) - bandQ*np.sin(dPhase)
    Q = bandI*np.sin(dPhase) + bandQ*np.cos(dPhase)
    SLCderamp = I + 1j*Q
    return SLCderamp


def getBandName(bandList,band):
    if len(band) < 2:
        out = [b for b in bandList if b.startswith(band)] # single string
    else:
        out = [b for b in bandList 
                if all(x in b for x in band)] # multiple sub-strings
    if out: # test if any matches
        return out[0]
    else:
        raise NameError('No band matching sub-strings:'+str(band))


def readSLC(file,boundingBox,method = 'coreg',deramp = False):
    """Read SLC data from SNAP file as np.array(range,azimuth) complex64
       and perform deramping & demodulation if requested
    
    input: .dim full file path
            boundingBox = ((minAz,maxAz),(minR,maxR))
            method = 'coreg' / 'raw'
    """    
    
    product = ProductIO.readProduct(file)
    bandNames = list(product.getBandNames())
    
    minAz = int(boundingBox[0][0])
    maxAz = int(boundingBox[0][1])
    minR = int(boundingBox[1][0])
    maxR = int(boundingBox[1][1])
    sizeAz = maxAz-minAz+1
    sizeR = maxR-minR+1

    # Initialise 
    array = np.zeros((sizeAz,sizeR),dtype=np.float32)

    if method == 'coreg':
        bandName = getBandName(bandNames,['i_','slv'])
        band = product.getBand(bandName)
        bandI = band.readPixels(minR, minAz, sizeR, sizeAz, array.copy())
        bandName = getBandName(bandNames,['q_','slv'])
        band = product.getBand(bandName)
        bandQ = band.readPixels(minR, minAz, sizeR, sizeAz, array.copy())

        if deramp:
            bandName = getBandName(bandNames,'derampDemodPhase_slv')
            band = product.getBand(bandName)
            dPhase = band.readPixels(minR, minAz, sizeR, sizeAz, array.copy())
            I = bandI*np.cos(dPhase) - bandQ*np.sin(dPhase)
            Q = bandI*np.sin(dPhase) + bandQ*np.cos(dPhase)
            SLC = I + 1j*Q
        else:
            SLC = bandI + 1j*bandQ   
    elif method == 'raw':
        bandName = getBandName(bandNames,'i_')
        band = product.getBand(bandName)
        bandI = band.readPixels(minR, minAz, sizeR, sizeAz, array.copy())
        bandName = getBandName(bandNames,'q_')
        band = product.getBand(bandName)
        bandQ = band.readPixels(minR, minAz, sizeR, sizeAz, array.copy())
        
        if deramp:
            dPhase = getDerampDemodPhase(file,boundingBox)
            I = bandI*np.cos(dPhase) - bandQ*np.sin(dPhase)
            Q = bandI*np.sin(dPhase) + bandQ*np.cos(dPhase)
            SLC = I + 1j*Q    
        else:
            SLC = bandI + 1j*bandQ

    product.dispose()    
    return SLC


def readDerampPhase(file,boundingBox):

    product = ProductIO.readProduct(file)
    bandNames = list(product.getBandNames())    
    minAz = int(boundingBox[0][0])
    maxAz = int(boundingBox[0][1])
    minR = int(boundingBox[1][0])
    maxR = int(boundingBox[1][1])
    sizeAz = maxAz-minAz+1
    sizeR = maxR-minR+1
    array = np.zeros((sizeAz,sizeR),dtype=np.float32)    
    bandName = getBandName(bandNames,'derampDemodPhase_slv')
    band = product.getBand(bandName)
    dPhase = band.readPixels(minR, minAz, sizeR, sizeAz, array.copy())
    product.dispose()
    return dPhase


def readIFG(file,boundingBox):
    
    # get ifg file corresponding to stack:
    file = file.replace('stack','itfg')
    
    product = ProductIO.readProduct(file)
    bandNames = list(product.getBandNames())

    minAz = int(boundingBox[0][0])
    maxAz = int(boundingBox[0][1])
    minR = int(boundingBox[1][0])
    maxR = int(boundingBox[1][1])
    sizeAz = maxAz-minAz+1
    sizeR = maxR-minR+1

    # Initialise 
    array = np.zeros((sizeAz,sizeR),dtype=np.float64)
    
    bandName = getBandName(bandNames,'i_ifg')
    band = product.getBand(bandName)
    bandI = band.readPixels(minR, minAz, sizeR, sizeAz, array.copy())
    bandName = getBandName(bandNames,'q_ifg')
    band = product.getBand(bandName)
    bandQ = band.readPixels(minR, minAz, sizeR, sizeAz, array.copy())
    IFG = bandI + 1j*bandQ
    
    product.dispose()
    return IFG


def readHGT(file,boundingBox):
    
    # get ifg file corresponding to stack:
    file = file.replace('stack','itfg')
    
    product = ProductIO.readProduct(file)
    bandNames = list(product.getBandNames())

    minAz = int(boundingBox[0][0])
    maxAz = int(boundingBox[0][1])
    minR = int(boundingBox[1][0])
    maxR = int(boundingBox[1][1])
    sizeAz = maxAz-minAz+1
    sizeR = maxR-minR+1

    # Initialise 
    array = np.zeros((sizeAz,sizeR),dtype=np.float64)
    
    if 'elevation' in bandNames:
        band = product.getBand('elevation')
        hgt = band.readPixels(minR, minAz, sizeR, sizeAz, array.copy())
        return hgt
    else:
        raise NameError('No elevation band in specified product.')
    product.dispose()


def bistaticCorrection(tRange,metadata):
    """Return a Sentinel-1 azimuth bistatic correction in [s]
    
    This fixes the S-1 IPF mid-IW2 "bulk" correction which is wrong and instead
    applies correct bistatic azimuth timing fix for given range.
    It is a correction = should be added to measured azimuth time!
    
    input: range time in [s] and metadata dict.
    """         

    rangeMid = round(metadata['nRange']/2)
    tRangeMid = metadata['range0time'] + (rangeMid-1)/metadata['RSR'];
    
    # fix to IW2 mid-range! (+ exclude null samples):
    if metadata['swath'] == 1:
        tRangeMid += (metadata['nRange']-350)/metadata['RSR']
    elif metadata['swath'] == 3:
        tRangeMid -= (metadata['nRange']-1680)/metadata['RSR']


    
    #bistaticAz = (tRange - tRangeMid)/2
    # first put back the IPF mid correction and then use range specific one:
    bistaticAzOLD = tRangeMid/2 + tRange/2 - metadata['range0time'] # The range0time is equal to rank*PRI only when the entire image is loaded
    
    # correctly should be (Balss, 2018):
    bistaticAz = tRangeMid + tRange - metadata['PRI']*metadata['rank']
    # PRI and rank in origMetadata/generalAnotation/downlink/

    # print('@BISTATIC VALIDATION --------------- ') # DEBUG ok
    # print(f'acqDate = {metadata["acqDate"]}')
    # print(f'rangeMid = {rangeMid}')
    # print(f'range0time = {metadata["range0time"]}')
    # print(f'RSR = {metadata["RSR"]}')
    # print(f'tRangeMid = {tRangeMid}')
    # print(f'tRange = {tRange}')
    # print(f'bistaticAz = {bistaticAz}')
    # print(f'rank*PRI = {metadata["PRI"]*metadata["rank"]}')
    # azimuthFactor = metadata['PRF']*metadata['azimuthSpacing']
    # print(f'bistaticAz [m] = {azimuthFactor*bistaticAz}')
    # print('------------------------------------ ')

    return bistaticAz


def dopplerRgCorrection(tRange,tAzimuth,metadata):
    
    #txPulseRampRate
    chirpRate = metadata['chirpRate']
    #Kr = 7.792817275120481e11    
    #B = 56.5e6 # S1 range bandwidth
    #pulseLength = 5.339056426703426e-05
    #Kr = B/pulseLength
    
    # get azimuth fm rate for all bursts:
    # mid-burst azimuth time
    if metadata['nBursts'] > 1:
        burstAz = metadata['burstInfo']['firstAzTime']
        idx = np.nonzero(np.mean(tAzimuth) > burstAz)[0][-1]
        tmp1 = metadata['burstInfo']['nAzBurst']/2/metadata['PRF'] # burst +-time span (Eta)
        midAzimuthTime = burstAz[idx] + tmp1 # burst +-time span (Eta)
    else:
        tmp1 = metadata['nAzimuth']/2/metadata['PRF'] # burst +-time span (Eta)
        midAzimuthTime = metadata['azimuth0time'] + tmp1 # burst +-time span (Eta)
        
    # velocity:
    dummy,vel = geoUtils.orbitVal(metadata['orbitFit'],midAzimuthTime)
    vel = np.linalg.norm(vel) # absolute value
    # angular steering rate:
    krot = 2*vel*metadata['steeringRate']/metadata['wavelength']
    
    # compute range dependent DopplerRate:
    azFmRateArray = metadata['azFmRateArray']
    tBurst = [np.mod(tmp['t']/3600/24,1)*3600*24 for tmp in azFmRateArray]
    
    # find appropriate burst index: 
    burstIdx = np.argmin(np.abs(tBurst-np.mean(tAzimuth)))
   
    dt = tRange*2 - azFmRateArray[burstIdx]['t0']
    rangeDependDopplerRate = (azFmRateArray[burstIdx]['c0'] + 
                              azFmRateArray[burstIdx]['c1']*dt + 
                              azFmRateArray[burstIdx]['c2']*dt*dt)
    
    # DopperRate [Hz/s]:
    dopplerRate = rangeDependDopplerRate*krot/(rangeDependDopplerRate - krot)
    
    #% get Doppler centroid:
    dcPolyArray = metadata['dcPolyArray']
    dt = tRange*2 - dcPolyArray[burstIdx]['t0']
    dopplerCentroid = (dcPolyArray[burstIdx]['c0'] + 
                       dcPolyArray[burstIdx]['c1']*dt + 
                       dcPolyArray[burstIdx]['c2']*dt*dt)
    # print(f'@@ DOPPLER validation ------------ on {str(metadata["acqDate"])[0:10]}') # DEBUG ok
    # print(f'dopplerCentroid = {dopplerCentroid}')
    # print(f'dopplerRate = {dopplerRate}')
    # print(f'tAzimuth = {tAzimuth}')
    # print(f'tmp1 = {tmp1}')
    # print(f'azimuth0time = {metadata["azimuth0time"]}')
    # print(f'nAzimuth = {metadata["nAzimuth"]}')
    # print(f'midAzimuthTime = {midAzimuthTime}')
    # print(f'chirpRate = {chirpRate}')
    # print(f'---------------------------------------------------')
    
    return (dopplerCentroid+dopplerRate*(tAzimuth-midAzimuthTime))/chirpRate/2

def FMmismatchCorrection(xyz,tRange,tAzimuth,metadata):
    
    # mid-burst azimuth time
    if metadata['nBursts'] > 1:
        burstAz = metadata['burstInfo']['firstAzTime']
        idx = np.nonzero(np.mean(tAzimuth) > burstAz)[0][-1]
        tmp1 = metadata['burstInfo']['nAzBurst']/2/metadata['PRF'] # burst +-time span (Eta)
        midAzimuthTime = burstAz[idx] + tmp1 # burst +-time span (Eta)
    else:
        tmp1 = metadata['nAzimuth']/2/metadata['PRF'] # burst +-time span (Eta)
        midAzimuthTime = metadata['azimuth0time'] + tmp1 # burst +-time span (Eta)
        
    # velocity:
    dummy,vel = geoUtils.orbitVal(metadata['orbitFit'],midAzimuthTime)
    vel = np.linalg.norm(vel) # absolute value
    # angular steering rate:
    krot = 2*vel*metadata['steeringRate']/metadata['wavelength']
    
    # compute range dependent DopplerRate:
    azFmRateArray = metadata['azFmRateArray']
    tBurst = [np.mod(tmp['t']/3600/24,1)*3600*24 for tmp in azFmRateArray]
    # find appropriate burst index: 
    burstIdx = np.argmin(np.abs(tBurst-np.mean(tAzimuth)))
   
    dt = tRange*2 - azFmRateArray[burstIdx]['t0']
    rangeDependDopplerRate = (azFmRateArray[burstIdx]['c0'] + 
                              azFmRateArray[burstIdx]['c1']*dt + 
                              azFmRateArray[burstIdx]['c2']*dt*dt)
    # DopperRate [Hz/s]:
    dopplerRate = rangeDependDopplerRate*krot/(rangeDependDopplerRate - krot)
    
    #% get Doppler centroid:
    dcPolyArray = metadata['dcPolyArray']
    dt = tRange*2 - dcPolyArray[burstIdx]['t0']
    dopplerCentroid = (dcPolyArray[burstIdx]['c0'] + 
                       dcPolyArray[burstIdx]['c1']*dt + 
                       dcPolyArray[burstIdx]['c2']*dt*dt)
    
    #referenceTIme
    midRangeTime = metadata['range0time']*2+metadata['nRange']/2/metadata['RSR']    
    dt0 = midRangeTime - dcPolyArray[burstIdx]['t0']
    dc0 = (dcPolyArray[burstIdx]['c0'] + dcPolyArray[burstIdx]['c1']*dt0 + 
       dcPolyArray[burstIdx]['c2']*dt0*dt0)
    dt0 = midRangeTime - azFmRateArray[burstIdx]['t0']
    fm0 = (azFmRateArray[burstIdx]['c0'] + azFmRateArray[burstIdx]['c1']*dt0 + 
       azFmRateArray[burstIdx]['c2']*dt0*dt0)
    # finally reference time:
    refTime = dopplerCentroid/rangeDependDopplerRate- dc0/fm0
    
    fDC = dopplerCentroid+dopplerRate*(tAzimuth-midAzimuthTime-refTime)
    
    # get true (geo) FM rate:
    satVec = radarUtils.xyz2satvec(xyz,metadata)  
    k_geo = -2/(metadata['wavelength']*np.linalg.norm(satVec[:3]-xyz))*(
        (satVec[:3]-xyz)@satVec[-3:] + satVec[3:-3]@satVec[3:-3])
    
    # calculate correction:
    return fDC*(1/(-rangeDependDopplerRate[0])-(1/(-k_geo)))
    

def inSLCframe(slcFile,wktAOI):
    
    import zipfile
    import re
    import xml.etree.ElementTree as ET
    from shapely import wkt,geometry
    
    ## debug
    #slcFile = '/data/data3/SNT1/2019/ASC102/S1A_IW_SLC__1SDV_20190718T162614_20190718T162641_028174_032EB2_343B.zip'
    #wktAOI = 'POLYGON ((21.34 48.60, 21.37 48.60, 21.37 48.69, 21.34 48.69, 21.34 48.60))'
    
    # open manifest.safe from SLC .zip file:
    try:
        zf = zipfile.ZipFile(slcFile, 'r')
    except:
        raise Exception("Zip file corrupt.")
    manifest = [n for n in zf.namelist() if re.search("manifest.safe",n)]
    f = zf.open(manifest[0])
    
    # parse gml frame from manifest.safe:
    root = ET.fromstring(f.read())
    frame = root[1].find("./metadataObject[@ID='measurementFrameSet']")
    gmlFrame = frame[0][0][0][0][0][0].text 
    
    # covert gml frame to shapely polygon:
    gmlFrame = [[float(i) for i in j.split(',')] for j in gmlFrame.split(' ')]
    frame = geometry.Polygon([i[::-1] for i in gmlFrame]) # change point order for shapely
    
    aoi = wkt.loads(wktAOI) # convert input wkt to shapely polygon
    
    return frame.contains(aoi)
