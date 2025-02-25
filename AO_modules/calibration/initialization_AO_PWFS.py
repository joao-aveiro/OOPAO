
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 16:41:21 2020

@author: cheritie
"""
# local modules 
from AO_modules.Telescope         import Telescope
from AO_modules.Source            import Source
from AO_modules.Atmosphere        import Atmosphere
from AO_modules.Pyramid           import Pyramid
from AO_modules.DeformableMirror  import DeformableMirror
from AO_modules.MisRegistration   import MisRegistration
from AO_modules.tools.set_paralleling_setup import set_paralleling_setup
import numpy as np
from AO_modules.calibration.compute_KL_modal_basis import compute_M2C
from AO_modules.calibration.ao_calibration import ao_calibration
from AO_modules.tools.tools import *

def run_initialization_AO_PWFS(param):
    
    #%% -----------------------     TELESCOPE   ----------------------------------
    
    # create the Telescope object
    tel = Telescope(resolution          = param['resolution'],\
                    diameter            = param['diameter'],\
                    samplingTime        = param['samplingTime'],\
                    centralObstruction  = param['centralObstruction'])
    
    #%% -----------------------     NGS   ----------------------------------
    # create the Source object
    ngs=Source(optBand   = param['opticalBand'],\
               magnitude = param['magnitude'])
    
    # combine the NGS to the telescope using '*' operator:
    ngs*tel
    
    
    
    #%% -----------------------     ATMOSPHERE   ----------------------------------

    # create the Atmosphere object
    atm=Atmosphere(telescope     = tel,\
                   r0            = param['r0'],\
                   L0            = param['L0'],\
                   windSpeed     = param['windSpeed'],\
                   fractionalR0  = param['fractionnalR0'],\
                   windDirection = param['windDirection'],\
                   altitude      = param['altitude'],\
                   param         = param)
    
    # initialize atmosphere
    atm.initializeAtmosphere(tel)
    
    # pairing the telescope with the atmosphere using '+' operator
    tel+atm
    # separating the tel and atmosphere using '-' operator
    tel-atm
    
    #%% -----------------------     DEFORMABLE MIRROR   ----------------------------------
    # mis-registrations object
    misReg = MisRegistration(param)
    # if no coordonates specified, create a cartesian dm
    dm=DeformableMirror(telescope    = tel,\
                        nSubap       = param['nSubaperture'],\
                        mechCoupling = param['mechanicalCoupling'],\
                        misReg       = misReg)
    
    
    #%% -----------------------     PYRAMID WFS   ----------------------------------
    
    # make sure tel and atm are separated to initialize the PWFS
    tel-atm
    try:
        print('Using a user-defined zero-padding for the PWFS:', param['zeroPadding'])
    except:
        param['zeroPadding']=0
        print('Using the default zero-padding for the PWFS')
    # create the Pyramid Object
    wfs = Pyramid(nSubap                = param['nSubaperture'],\
                  telescope             = tel,\
                  modulation            = param['modulation'],\
                  lightRatio            = param['lightThreshold'],\
                  pupilSeparationRatio  = param['pupilSeparationRatio'],\
                  calibModulation       = param['calibrationModulation'],\
                  psfCentering          = param['psfCentering'],\
                  edgePixel             = param['edgePixel'],\
                  postProcessing        = param['postProcessing'])    
    # set_paralleling_setup(wfs, ELT = False)
    # propagate the light through the WFS
    tel*wfs
    #%% -----------------------     Modal Basis   ----------------------------------

    try:
        M2C_KL = read_fits(param['pathInput']+param['modal_basis_name']+'.fits')
        print('Succesfully loaded KL modal basis')
    except:
        print('Computing KL modal basis ...')

        M2C_KL = compute_M2C(    telescope          = tel,\
                                 atmosphere         = atm,\
                                 deformableMirror   = dm,\
                                 param              = param,\
                                 nameFolder         = None,\
                                 nameFile           = None,\
                                 remove_piston      = True,\
                                 HHtName            = None,\
                                 baseName           = None ,\
                                 mem_available      = 8.1e9,\
                                 minimF             = False,\
                                 nmo                = 1000,\
                                 ortho_spm          = True,\
                                 SZ                 = np.int(2*tel.OPD.shape[0]),\
                                 nZer               = 3,\
                                 NDIVL              = 1)
        print('Done')
    

#%% -----------------------     Calibration   ----------------------------------

    ao_calib =  ao_calibration(param            = param,\
                               ngs              = ngs,\
                               tel              = tel,\
                               atm              = atm,\
                               dm               = dm,\
                               wfs              = wfs,\
                               nameFolderIntMat = None,\
                               nameIntMat       = None,\
                               nameFolderBasis  = None,\
                               nameBasis        = param['modal_basis_name'],\
                               nMeasurements    = 1,\
                               get_basis        = False)
    
   
        

    class emptyClass():
        pass
    # save output as sub-classes
    simulationObject = emptyClass()
    
    simulationObject.tel    = tel
    simulationObject.atm    = atm

    simulationObject.ngs    = ngs
    simulationObject.dm     = dm
    simulationObject.wfs    = wfs
    simulationObject.param  = param
    simulationObject.calib  = ao_calib.calib
    simulationObject.M2C    = ao_calib.M2C
    simulationObject.gOpt   = ao_calib.gOpt
    simulationObject.projector   = ao_calib.projector

    return simulationObject
