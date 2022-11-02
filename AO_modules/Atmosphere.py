# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 10:59:02 2020

@author: cheritie
"""
import inspect
import numpy             as np
from AO_modules.phaseStats import makeCovarianceMatrix,ft_phase_screen
from AO_modules.tools.tools import emptyClass,translationImageMatrix,globalTransformation, createFolder, pol2cart, cart2pol
from AO_modules.tools.interpolateGeometricalTransformation import interpolate_cube
from AO_modules.tools.displayTools import getColorOrder

import time
import json
import jsonpickle
import matplotlib.pyplot as plt
from numpy.random import RandomState
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CLASS INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
class Atmosphere:
    def __init__(self,telescope,r0,L0,windSpeed,fractionalR0,windDirection,altitude,mode=2, param = None, asterism = None):
        """
        An Atmosphere object is made of one or several layer of turbulence that follow the Van Karmann statistics. Each layer is considered to be independant to the other ones and has its own properties (direction, speed, etc.)
        
        The Atmosphere object can be defined for a single Source object (default) or multi Source Object. The Source coordinates allow to span different areas in the field (defined as well by the tel.fov).
        
        If the source type is an LGS the cone effect is considered using an interpolation. NGS and LGS can be combined together in the asterism object. 
        
        The convention chosen is that all the wavelength-dependant atmosphere parameters are expressed at 500 nm. 

        
        ************************** REQUIRED PARAMETERS **************************
        It requires the following parameters: 
            
        _ telescope             : the telescope object to which the Atmosphere is associated. This object carries the phase, flux, pupil information and sampling time as well as the type of source (NGS/LGS, source/asterism)
        _ r0                    : the Fried parameter in m 
        _ L0                    : Outer scale parameter
        _ windSpeed             : List of wind-speed for each layer in [m/s]
        _ fractionalR0          : Cn2 profile of the turbulence. This should be a list of values for each layer
        _ windDirection         : List of wind-direction for each layer in [deg]
        _ altitude              : List of altitude for each layer in [m]
        
        ************************** OPTIONAL PARAMETERS **************************
        
        _ mode                  : Different ways to generate the spectrum of the atmosphere
        _ param                 : parameter file of the system. Once computed, the covariance matrices are saved in the calibration data folder and loaded instead of re-computed evry time.
        _ asterism              : if the system contains multiple source, an astrism should be input to the atmosphere object
        
        ************************** COUPLING A TELESCOPE AND AN ATMOSPHERE OBJECT **************************
        A Telescope object "tel" can be coupled to an Atmosphere object "atm" using: 
            _ tel + atm
        This means that a bridge is created between atm and tel: everytime that atm.OPD is updated, the tel.OPD property is automatically set to atm.OPD to reproduce the effect of the turbulence.
        
        A Telescope object "tel" can be separated of an Atmosphere object "atm" using: 
            _ tel - atm
        This corresponds to a diffraction limited case (no turbulence)

    
        ************************** PROPERTIES **************************
        
        The main properties of the Atmosphere object are listed here: 
    
        _ atm.OPD : Optical Path Difference in [m] truncated by the telescope pupil. If the atmosphere has multiple sources, the OPD is a list of OPD for each source
        _ atm.OPD_no_pupil : Optical Path Difference in [m]. If the atmosphere has multiple sources, the OPD is a list of OPD for each source
        _ atm.r0                   
        _ atm.L0
        _ atm.nLayer                                : number of turbulence layers
        _ atm.seeingArcsec                          : seeing in arcsec at 500 nm
        _ atm.layer_X                               : access the child object corresponding to the layer X where X starts at 0   
        
        the following properties can be updated on the fly:
            _ atm.r0            
            _ atm.windSpeed      
            _ atm.windDirection 
        ************************** FUNCTIONS **************************

        _ atm.update()                              : update the OPD of the atmosphere for each layer according to the time step defined by tel.samplingTime                    
        _ atm.generateNewPhaseScreen(seed)          : generate a new phase screen for the atmosphere OPD
        _ atm.print_atm_at_wavelength(wavelength)   : prompt seeing and r0 at specified wavelength
        _ atm.print_atm()                           : prompt the main properties of the atm object
        _ display_atm_layers(layer_index)           : imshow the OPD of each layer with the intersection beam for each source
        
        """        
        self.hasNotBeenInitialized  = True
        self.r0_def                 = 0.15                # Fried Parameter in m 
        self.r0                     = r0                # Fried Parameter in m 
        self.fractionalR0           = fractionalR0      # Cn2 square profile
        self.L0                     = L0                # Outer Scale in m
        self.altitude               = altitude          # altitude of the layers
        self.nLayer                 = len(fractionalR0)     # number of layer
        self.windSpeed              = windSpeed         # wind speed of the layers in m/s
        self.windDirection          = windDirection     # wind direction in degrees
        self.tag                    = 'atmosphere'      # Tag of the object
        self.nExtra                 = 2                 # number of extra pixel to generate the phase screens
        self.wavelength             = 500*1e-9          # Wavelengt used to define the properties of the atmosphere
        self.tel                    = telescope         # associated telescope object
        self.mode                   = mode              # DEBUG -> first phase screen generation mode
        self.seeingArcsec           = 206265*(self.wavelength/self.r0)
        self.asterism               = asterism          # case when multiple sources are considered (LGS and NGS)
        self.param                  = param
        if self.asterism is not None:
            self.oversampling_factor    = np.max((np.asarray(self.asterism.coordinates)[:,0]/(self.tel.resolution/2)))
        else:
            self.oversampling_factor = self.tel.src.coordinates[0]/(self.tel.resolution/2)
        
        
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ATM INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    def initializeAtmosphere(self,telescope):
        phase_support = self.initialize_phase_support()
        if self.hasNotBeenInitialized:
            self.initial_r0 = self.r0
            for i_layer in range(self.nLayer):       
                # create the layer
                print('Creation of layer' + str(i_layer+1) + '/' + str(self.nLayer) + ' ...' )
                tmpLayer=self.buildLayer(telescope,self.r0_def,self.L0,i_layer = i_layer)
                setattr(self,'layer_'+str(i_layer+1),tmpLayer)
                
                phase_support = self.fill_phase_support(tmpLayer,phase_support,i_layer)
                tmpLayer.phase_support = phase_support
                tmpLayer.phase *= self.wavelength/2/np.pi

        else:
            print('Re-setting the atmosphere to its initial state...' )
            self.r0 = self.initial_r0
            for i_layer in range(self.nLayer):       
                print('Updating layer' + str(i_layer+1) + '/' + str(self.nLayer) + ' ...' )
                tmpLayer = getattr(self,'layer_'+str(i_layer+1))
                tmpLayer.phase          = tmpLayer.initialPhase/self.wavelength*2*np.pi
                tmpLayer.randomState    = RandomState(42+i_layer*1000)
               
                Z = tmpLayer.phase[tmpLayer.innerMask[1:-1,1:-1]!=0]
                X = np.matmul(tmpLayer.A,Z) + np.matmul(tmpLayer.B,tmpLayer.randomState.normal( size=tmpLayer.B.shape[1]))
                
                tmpLayer.mapShift[tmpLayer.outerMask!=0] = X
                tmpLayer.mapShift[tmpLayer.outerMask==0] = np.reshape(tmpLayer.phase,tmpLayer.resolution*tmpLayer.resolution)
                
                tmpLayer.notDoneOnce = True                

                setattr(self,'layer_'+str(i_layer+1),tmpLayer) 
                
                phase_support = self.fill_phase_support(tmpLayer,phase_support,i_layer)

                # wavelenfth scaling
                tmpLayer.phase *= self.wavelength/2/np.pi
        self.generateNewPhaseScreen(seed=0)
                
        self.hasNotBeenInitialized  = False        
        # save the resulting phase screen in OPD  
        self.set_OPD(phase_support)
        self.print_atm()
            
    def buildLayer(self,telescope,r0,L0,i_layer,compute_turbulence = True):
        """
            Generation of phase screens using the method introduced in Assemat et al (2006)
        """

    
        # initialize layer object
        layer               = emptyClass()
        # create a random state to allow reproductible sequences of phase screens
        layer.randomState   = RandomState(42+i_layer*1000)
        
        # gather properties of the atmosphere
        layer.altitude      = self.altitude[i_layer]       
        layer.windSpeed     = self.windSpeed[i_layer]
        layer.direction     = self.windDirection[i_layer]
        
        # compute the X and Y wind speed
        layer.vY            = layer.windSpeed*np.cos(np.deg2rad(layer.direction))
        layer.vX            = layer.windSpeed*np.sin(np.deg2rad(layer.direction))      
        
        # Diameter and resolution of the layer including the Field Of View and the number of extra pixels
        layer.D             = self.tel.D+2*np.tan(self.tel.fov/2)*layer.altitude*self.oversampling_factor
        layer.resolution    = int(np.ceil((self.tel.resolution/self.tel.D)*layer.D))
        
        layer.D_fov             = self.tel.D+2*np.tan(self.tel.fov/2)*layer.altitude
        layer.resolution_fov    = int(np.ceil((self.tel.resolution/self.tel.D)*layer.D))
        
        layer.center = layer.resolution//2
        
        if self.asterism is None:
            [x_z,y_z] = pol2cart(self.tel.src.coordinates[0]*(layer.D_fov-self.tel.D)/self.tel.D,np.deg2rad(self.tel.src.coordinates[1]))
     
            center_x = int(y_z)+layer.resolution//2
            center_y = int(x_z)+layer.resolution//2
        
            layer.pupil_footprint = np.zeros([layer.resolution,layer.resolution])
            layer.pupil_footprint[center_x-self.tel.resolution//2:center_x+self.tel.resolution//2,center_y-self.tel.resolution//2:center_y+self.tel.resolution//2 ] = 1
        else:
            layer.pupil_footprint= []
            for i in range(self.asterism.n_source):
                 [x_z,y_z] = pol2cart(self.asterism.coordinates[i][0]*(layer.D_fov-self.tel.D)/self.tel.D,np.deg2rad(self.asterism.coordinates[i][1]))
     
                 center_x = int(y_z)+layer.resolution//2
                 center_y = int(x_z)+layer.resolution//2
            
                 pupil_footprint = np.zeros([layer.resolution,layer.resolution])
                 pupil_footprint[center_x-self.tel.resolution//2:center_x+self.tel.resolution//2,center_y-self.tel.resolution//2:center_y+self.tel.resolution//2 ] = 1
                 layer.pupil_footprint.append(pupil_footprint)   
        # layer pixel size
        layer.d0            = layer.D/layer.resolution
        
        # number of pixel for the phase screens computation
        layer.nExtra        = self.nExtra
        layer.nPixel        = int(1+np.round(layer.D/layer.d0))
        if compute_turbulence:
            print('-> Computing the initial phase screen...')  
            a=time.time()
            if self.mode ==1:
                import aotools as ao
                layer.phaseScreen   = ao.turbulence.infinitephasescreen.PhaseScreenVonKarman(layer.resolution,layer.D/(layer.resolution),r0,L0,random_seed=i_layer)
                layer.phase         = layer.phaseScreen.scrn
            else:
                if self.mode == 2:
                    from AO_modules.phaseStats import ft_sh_phase_screen
    
                    layer.phase         = ft_sh_phase_screen(self,layer.resolution,layer.D/layer.resolution,seed=i_layer)                
                else: 
                    layer.phase         = ft_phase_screen(self,layer.resolution,layer.D/layer.resolution,seed=i_layer)
                    
            layer.initialPhase = layer.phase
            layer.seed = i_layer
            b=time.time()
            print('initial phase screen : ' +str(b-a) +' s')
            
            # Outer ring of pixel for the phase screens update 
            layer.outerMask             = np.ones([layer.resolution+layer.nExtra,layer.resolution+layer.nExtra])
            layer.outerMask[1:-1,1:-1]  = 0
            
            # inner pixels that contains the phase screens
            layer.innerMask             = np.ones([layer.resolution+layer.nExtra,layer.resolution+layer.nExtra])
            layer.innerMask -= layer.outerMask
            layer.innerMask[1+layer.nExtra:-1-layer.nExtra,1+layer.nExtra:-1-layer.nExtra] = 0
            
            l = np.linspace(0,layer.resolution+1,layer.resolution+2) * layer.D/(layer.resolution-1)
            u,v = np.meshgrid(l,l)
            
            layer.innerZ = u[layer.innerMask!=0] + 1j*v[layer.innerMask!=0]
            layer.outerZ = u[layer.outerMask!=0] + 1j*v[layer.outerMask!=0]
            
            self.get_covariance_matrices(layer)                            
            
            layer.A         = np.matmul(self.ZXt_r0.T,self.ZZt_inv_r0)    
            BBt             = self.XXt_r0 -  np.matmul(layer.A,self.ZXt_r0)
            layer.B         = np.linalg.cholesky(BBt)
            layer.mapShift  = np.zeros([layer.nPixel+1,layer.nPixel+1])        
            Z               = layer.phase[layer.innerMask[1:-1,1:-1]!=0]
            X               = np.matmul(layer.A,Z) + np.matmul(layer.B,layer.randomState.normal(size=layer.B.shape[1]))
            
            layer.mapShift[layer.outerMask!=0] = X
            layer.mapShift[layer.outerMask==0] = np.reshape(layer.phase,layer.resolution*layer.resolution)
            layer.notDoneOnce                  = True

        print('Done!')
        
    
        return layer
    
    
    def add_row(self,layer,stepInPixel):
        shiftMatrix                         = translationImageMatrix(layer.mapShift,[stepInPixel[0],stepInPixel[1]]) #units are in pixel of the M1            
        tmp                                 = globalTransformation(layer.mapShift,shiftMatrix)
        onePixelShiftedPhaseScreen          = tmp[1:-1,1:-1]        
        Z                                   = onePixelShiftedPhaseScreen[layer.innerMask[1:-1,1:-1]!=0]
        X                                   = layer.A@Z + layer.B@layer.randomState.normal(size=layer.B.shape[1])
        layer.mapShift[layer.outerMask!=0]  = X
        layer.mapShift[layer.outerMask==0]  = np.reshape(onePixelShiftedPhaseScreen,layer.resolution*layer.resolution)
        return onePixelShiftedPhaseScreen

    def updateLayer(self,layer):
        self.ps_loop    = layer.D / (layer.resolution)
        ps_turb_x       = layer.vX*self.tel.samplingTime
        ps_turb_y       = layer.vY*self.tel.samplingTime
        
        if layer.vX==0 and layer.vY==0:
            layer.phase = layer.phase
            
        else:
            if layer.notDoneOnce:
                layer.notDoneOnce = False
                layer.ratio = np.zeros(2)
                layer.ratio[0] = ps_turb_x/self.ps_loop
                layer.ratio[1] = ps_turb_y/self.ps_loop
                layer.buff = np.zeros(2)
            tmpRatio = np.abs(layer.ratio)
            tmpRatio[np.isinf(tmpRatio)]=0
            nScreens = (tmpRatio)
            nScreens = nScreens.astype('int')
            
            stepInPixel =np.zeros(2)
            stepInSubPixel =np.zeros(2)
            
            for i in range(nScreens.min()):
                stepInPixel[0]=1
                stepInPixel[1]=1
                stepInPixel=stepInPixel*np.sign(layer.ratio)
                layer.phase = self.add_row(layer,stepInPixel)
                
            for j in range(nScreens.max()-nScreens.min()):   
                stepInPixel[0]=1
                stepInPixel[1]=1
                stepInPixel=stepInPixel*np.sign(layer.ratio)
                stepInPixel[np.where(nScreens==nScreens.min())]=0
                layer.phase = self.add_row(layer,stepInPixel)
            
            
            stepInSubPixel[0] =  (np.abs(layer.ratio[0])%1)*np.sign(layer.ratio[0])
            stepInSubPixel[1] =  (np.abs(layer.ratio[1])%1)*np.sign(layer.ratio[1])
            
            layer.buff += stepInSubPixel
            if np.abs(layer.buff[0])>=1 or np.abs(layer.buff[1])>=1:   
                stepInPixel[0] = 1*np.sign(layer.buff[0])
                stepInPixel[1] = 1*np.sign(layer.buff[1])
                stepInPixel[np.where(np.abs(layer.buff)<1)]=0    
                
                layer.phase = self.add_row(layer,stepInPixel)
    
            layer.buff[0]   =  (np.abs(layer.buff[0])%1)*np.sign(layer.buff[0])
            layer.buff[1]   =  (np.abs(layer.buff[1])%1)*np.sign(layer.buff[1])
                
            shiftMatrix     = translationImageMatrix(layer.mapShift,[layer.buff[0],layer.buff[1]]) #units are in pixel of the M1            
            layer.phase     = globalTransformation(layer.mapShift,shiftMatrix)[1:-1,1:-1]

    def update(self):
        phase_support = self.initialize_phase_support()
        for i_layer in range(self.nLayer):
            tmpLayer=getattr(self,'layer_'+str(i_layer+1))
            self.updateLayer(tmpLayer)
            phase_support = self.fill_phase_support(tmpLayer,phase_support,i_layer)
        self.set_OPD(phase_support)
        if self.tel.isPaired:
            self*self.tel
            
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ATM METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    def initialize_phase_support(self):
        if self.asterism is None:
            phase_support = np.zeros([self.tel.resolution,self.tel.resolution])
        else:
            phase_support = []
            for i in range(self.asterism.n_source):
                phase_support.append(np.zeros([self.tel.resolution,self.tel.resolution]))
        return phase_support
    def fill_phase_support(self,tmpLayer,phase_support,i_layer):
        if self.asterism is None:
            phase_support+= np.reshape(tmpLayer.phase[np.where(tmpLayer.pupil_footprint==1)],[self.tel.resolution,self.tel.resolution])* np.sqrt(self.fractionalR0[i_layer])
            # wavelenfth scaling
            tmpLayer.phase *= self.wavelength/2/np.pi
        else:
            for i in range(self.asterism.n_source):
                if self.asterism.src[i].type == 'LGS':
                    sub_im = np.reshape(tmpLayer.phase[np.where(tmpLayer.pupil_footprint[i]==1)],[self.tel.resolution,self.tel.resolution])
                    alpha_cone = np.arctan(self.tel.D/2/self.asterism.altitude[i])
                    h = self.asterism.altitude[i]-tmpLayer.altitude
                    if np.isinf(h):
                        r =self.tel.D/2
                    else:
                        r = h*np.tan(alpha_cone)
                    ratio = self.tel.D/r/2
                    cube_in = np.atleast_3d(sub_im).T
                    
                    pixel_size_in   = tmpLayer.D/tmpLayer.resolution
                    pixel_size_out  = pixel_size_in/ratio
                    resolution_out  = self.tel.resolution

                    phase_support[i]+= np.squeeze(interpolate_cube(cube_in, pixel_size_in, pixel_size_out, resolution_out)).T* np.sqrt(self.fractionalR0[i_layer])

                else:
                    phase_support[i]+= np.reshape(tmpLayer.phase[np.where(tmpLayer.pupil_footprint[i]==1)],[self.tel.resolution,self.tel.resolution])* np.sqrt(self.fractionalR0[i_layer])
        return phase_support
    
    def set_OPD(self,phase_support):
        if self.asterism is None:
            self.OPD_no_pupil   = phase_support*self.wavelength/2/np.pi
            self.OPD            = self.OPD_no_pupil*self.tel.pupil
        else:
            self.OPD =[]
            self.OPD_no_pupil = []
            for i in range(self.asterism.n_source):
                self.OPD_no_pupil.append(phase_support[i]*self.wavelength/2/np.pi)
                self.OPD.append(self.OPD_no_pupil[i]*self.tel.pupil)
        return
    def get_covariance_matrices(self,layer):
        # Compute the covariance matrices
        compute_covariance_matrices = True

        if self.tel.fov == 0: 
            try:
                c=time.time()        
                self.ZZt_r0 = self.ZZt_r0
                d=time.time()
                print('ZZt.. : ' +str(d-c) +' s')
                self.ZXt_r0 = self.ZXt_r0
                e=time.time()
                print('ZXt.. : ' +str(e-d) +' s')
                self.XXt_r0 = self.XXt_r0
                f=time.time()
                print('XXt.. : ' +str(f-e) +' s')
                self.ZZt_inv_r0 = self.ZZt_inv_r0
    
                print('SCAO system considered: covariance matrices were already computed!')
                compute_covariance_matrices = False
            except:                        
                compute_covariance_matrices = True
        if compute_covariance_matrices:
            c=time.time()        
            self.ZZt = makeCovarianceMatrix(layer.innerZ,layer.innerZ,self)
            
            if self.param is None:
                self.ZZt_inv = np.linalg.pinv(self.ZZt)
            else:
                try:
                    print('Loading pre-computed data...')            
                    name_data       = 'ZZt_inv_spider_L0_'+str(self.L0)+'_m_r0_'+str(self.r0_def)+'_shape_'+str(self.ZZt.shape[0])+'x'+str(self.ZZt.shape[1])+'.json'
                    location_data   = self.param['pathInput'] + self.param['name'] + '/sk_v/'
                    try:
                        with open(location_data+name_data ) as f:
                            C = json.load(f)
                        data_loaded = jsonpickle.decode(C)               
                    except:
                        createFolder(location_data)
                        with open(location_data+name_data ) as f:
                            C = json.load(f)
                        data_loaded = jsonpickle.decode(C)                    
                    self.ZZt_inv = data_loaded['ZZt_inv']
                    
                except: 
                    print('Something went wrong.. re-computing ZZt_inv ...')
                    name_data       = 'ZZt_inv_spider_L0_'+str(self.L0)+'_m_r0_'+str(self.r0_def)+'_shape_'+str(self.ZZt.shape[0])+'x'+str(self.ZZt.shape[1])+'.json'
                    location_data   = self.param['pathInput'] + self.param['name'] + '/sk_v/'
                    createFolder(location_data)
                    
                    self.ZZt_inv = np.linalg.pinv(self.ZZt)
                
                    print('saving for future...')
                    data = dict()
                    data['pupil'] = self.tel.pupil
                    data['ZZt_inv'] = self.ZZt_inv
                            
                    data_encoded  = jsonpickle.encode(data)
                    with open(location_data+name_data, 'w') as f:
                        json.dump(data_encoded, f)
            d=time.time()
            print('ZZt.. : ' +str(d-c) +' s')
            self.ZXt = makeCovarianceMatrix(layer.innerZ,layer.outerZ,self)
            e=time.time()
            print('ZXt.. : ' +str(e-d) +' s')
            self.XXt = makeCovarianceMatrix(layer.outerZ,layer.outerZ,self)
            f=time.time()
            print('XXt.. : ' +str(f-e) +' s')
            
            self.ZZt_r0     = self.ZZt*(self.r0_def/self.r0)**(5/3)
            self.ZXt_r0     = self.ZXt*(self.r0_def/self.r0)**(5/3)
            self.XXt_r0     = self.XXt*(self.r0_def/self.r0)**(5/3)
            self.ZZt_inv_r0 = self.ZZt_inv/((self.r0_def/self.r0)**(5/3))
        return
        
    def generateNewPhaseScreen(self,seed = None):
        if seed is None:
            t = time.localtime()
            seed = t.tm_hour*3600 + t.tm_min*60 + t.tm_sec
        phase_support = self.initialize_phase_support()
        for i_layer in range(self.nLayer):
            tmpLayer=getattr(self,'layer_'+str(i_layer+1))
            
            if self.mode ==1:
                import aotools as ao
                phaseScreen   = ao.turbulence.infinitephasescreen.PhaseScreenVonKarman(tmpLayer.resolution,tmpLayer.D/(tmpLayer.resolution),self.r0,self.L0,random_seed=seed+i_layer)
                phase         = phaseScreen.scrn
            else:
                if self.mode == 2:
                    from AO_modules.phaseStats import ft_sh_phase_screen
                    phase         = ft_sh_phase_screen(self,tmpLayer.resolution,tmpLayer.D/tmpLayer.resolution,seed=seed+i_layer)                
                else: 
                    phase         = ft_phase_screen(self,tmpLayer.resolution,tmpLayer.D/tmpLayer.resolution,seed=seed+i_layer)
            
            tmpLayer.phase = phase
            tmpLayer.randomState    = RandomState(42+i_layer*1000)
            
            Z = tmpLayer.phase[tmpLayer.innerMask[1:-1,1:-1]!=0]
            X = np.matmul(tmpLayer.A,Z) + np.matmul(tmpLayer.B,tmpLayer.randomState.normal( size=tmpLayer.B.shape[1]))
            
            tmpLayer.mapShift[tmpLayer.outerMask!=0] = X
            tmpLayer.mapShift[tmpLayer.outerMask==0] = np.reshape(tmpLayer.phase,tmpLayer.resolution*tmpLayer.resolution)
            tmpLayer.notDoneOnce = True

            setattr(self,'layer_'+str(i_layer+1),tmpLayer )
            phase_support = self.fill_phase_support(tmpLayer,phase_support,i_layer)
        self.set_OPD(phase_support)
        if self.tel.isPaired:
            self*self.tel


    def print_atm_at_wavelength(self,wavelength):

        r0_wvl              = self.r0*((wavelength/self.wavelength)**(5/3))
        seeingArcsec_wvl    = 206265*(wavelength/r0_wvl)

        print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ATMOSPHERE AT '+str(wavelength)+' nm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        print('r0 \t\t'+str(r0_wvl) + ' \t [m]') 
        print('Seeing \t' + str(np.round(seeingArcsec_wvl,2)) + str('\t ["]'))
        print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        return
        
    def print_atm(self):
        print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ATMOSPHERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        print('{: ^18s}'.format('Layer') + '{: ^18s}'.format('Direction [deg]')+ '{: ^18s}'.format('Speed [m/s]')+ '{: ^18s}'.format('Altitude [m]')+ '{: ^18s}'.format('Cn2 [m-2/3]') )
        print('------------------------------------------------------------------------------------------')
        
        for i_layer in range(self.nLayer):
            print('{: ^18s}'.format(str(i_layer+1)) + '{: ^18s}'.format(str(self.windDirection[i_layer]))+ '{: ^18s}'.format(str(self.windSpeed[i_layer]))+ '{: ^18s}'.format(str(self.altitude[i_layer]))+ '{: ^18s}'.format(str(self.fractionalR0[i_layer]) ))
            print('------------------------------------------------------------------------------------------')
        print('******************************************************************************************')

        print('{: ^18s}'.format('r0') + '{: ^18s}'.format(str(self.r0)+' [m]' ))
        print('{: ^18s}'.format('L0') + '{: ^18s}'.format(str(self.L0)+' [m]' ))
        print('{: ^18s}'.format('Seeing (V)') + '{: ^18s}'.format(str(np.round(self.seeingArcsec,2))+' ["]'))
        print('{: ^18s}'.format('Frequency') + '{: ^18s}'.format(str(np.round(1/self.tel.samplingTime,2))+' [Hz]' ))
        print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        
    def __mul__(self,obj):
        obj.OPD=self.OPD
        obj.OPD_no_pupil=self.OPD_no_pupil
        obj.isPaired=True
        return obj
    
    
    def display_atm_layers(self,layer_index,fig_index = None):
        if layer_index is None:
            layer_index = list(np.arange(self.nLayer))
        if type(layer_index) is not list:
            raise TypeError(' layer_index should be a list') 
        col = getColorOrder() 
        if fig_index is None:
            
            fig_index = time.time_ns()
        plt.figure(fig_index)
        axis_list = []
        for i in range(len(layer_index)):
            axis_list.append(plt.subplot(1,len(layer_index),i+1))
                             
        for i_l,ax in enumerate(axis_list):
            
            tmpLayer = getattr(self, 'layer_'+str(layer_index[i_l]+1))
            ax.imshow(tmpLayer.phase,extent = [-tmpLayer.D/2,tmpLayer.D/2,-tmpLayer.D/2,tmpLayer.D/2])
            center = tmpLayer.D/2
            [x_tel,y_tel] = pol2cart(tmpLayer.D_fov/2, np.linspace(0,2*np.pi,100,endpoint=True))  
            if self.tel.src.tag =='asterism':
                for i_source in range(len(self.tel.src.src)):
                    
                    [x_c,y_c] = pol2cart(self.tel.D/2, np.linspace(0,2*np.pi,100,endpoint=True))
                    alpha_cone = np.arctan(self.tel.D/2/self.tel.src.src[i_source].altitude)
                    
                    h = self.tel.src.src[i_source].altitude-tmpLayer.altitude
                    if np.isinf(h):
                        r =self.tel.D/2
                    else:
                        r = h*np.tan(alpha_cone)
                    [x_cone,y_cone] = pol2cart(r, np.linspace(0,2*np.pi,100,endpoint=True))
                    
                    [x_z,y_z] = pol2cart(self.tel.src.src[i_source].coordinates[0]*(tmpLayer.D_fov-self.tel.D)/self.tel.resolution,np.deg2rad(self.tel.src.src[i_source].coordinates[1]))
                    center = 0
                    [x_c,y_c] = pol2cart(tmpLayer.D_fov/2, np.linspace(0,2*np.pi,100,endpoint=True))  
                
                    ax.plot(x_cone+x_z+center,y_cone+y_z+center,'-', color = col [i_source])        
                    ax.fill(x_cone+x_z+center,y_cone+y_z+center,y_z+center, alpha = 0.25, color = col[i_source])
            else:
                [x_c,y_c] = pol2cart(self.tel.D/2, np.linspace(0,2*np.pi,100,endpoint=True))
                alpha_cone = np.arctan(self.tel.D/2/self.tel.src.altitude)
                
                h = self.tel.src.altitude-tmpLayer.altitude
                if np.isinf(h):
                    r =self.tel.D/2
                else:
                    r = h*np.tan(alpha_cone)
                [x_cone,y_cone] = pol2cart(r, np.linspace(0,2*np.pi,100,endpoint=True))
                
                [x_z,y_z] = pol2cart(self.tel.src.coordinates[0]*(tmpLayer.D_fov-self.tel.D)/self.tel.resolution,np.deg2rad(self.tel.src.coordinates[1]))
            
                center = 0
                [x_c,y_c] = pol2cart(tmpLayer.D_fov/2, np.linspace(0,2*np.pi,100,endpoint=True))  
            
                ax.plot(x_cone+x_z+center,y_cone+y_z+center,'-')        
                ax.fill(x_cone+x_z+center,y_cone+y_z+center,y_z+center, alpha = 0.6)     
            
            ax.set_xlabel('[m]')
            ax.set_ylabel('[m]')
            ax.set_title('Altitude '+str(tmpLayer.altitude)+' m')
            ax.plot(x_tel+center,y_tel+center,'--',color = 'k')
            
 # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ATM PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    
    @property
    def r0(self):
         return self._r0
    
    @r0.setter
    def r0(self,val):
         self._r0 = val

         if self.hasNotBeenInitialized is False:
             print('Updating the Atmosphere covariance matrices...')             
             self.ZZt_r0 = self.ZZt*(self.r0_def/self.r0)**(5/3)
             self.ZXt_r0 = self.ZXt*(self.r0_def/self.r0)**(5/3)
             self.XXt_r0 = self.XXt*(self.r0_def/self.r0)**(5/3)
             self.ZZt_inv_r0 = self.ZZt_inv/((self.r0_def/self.r0)**(5/3))
             
             self.seeingArcsec           = 206265*(self.wavelength/val)             
             for i_layer in range(self.nLayer):
                    tmpLayer = getattr(self,'layer_'+str(i_layer+1))
                    BBt                = self.XXt_r0 -  np.matmul(tmpLayer.A,self.ZXt_r0)
                    tmpLayer.B         = np.linalg.cholesky(BBt)

    @property
    def L0(self):
         return self._L0
    
    @L0.setter
    def L0(self,val):
         self._L0 = val
         if self.hasNotBeenInitialized is False:
             print('Updating the Atmosphere covariance matrices...')

             self.hasNotBeenInitialized = True
             del self.ZZt
             del self.XXt
             del self.ZXt
             del self.ZZt_inv
             self.initializeAtmosphere(self.tel)
    @property
    def windSpeed(self):
         return self._windSpeed
    
    @windSpeed.setter
    def windSpeed(self,val):
        self._windSpeed = val

        if self.hasNotBeenInitialized is False:
            if len(val)!= self.nLayer:
                print('Error! Wrong value for the wind-speed! Make sure that you inpute a wind-speed for each layer')
            else:
                print('Updating the wing speed...')
                for i_layer in range(self.nLayer):
                    tmpLayer = getattr(self,'layer_'+str(i_layer+1))
                    # tmpLayer.notDoneOnce = True

                    tmpLayer.windSpeed = val[i_layer]
                    tmpLayer.vY            = tmpLayer.windSpeed*np.cos(np.deg2rad(tmpLayer.direction))                    
                    tmpLayer.vX            = tmpLayer.windSpeed*np.sin(np.deg2rad(tmpLayer.direction))
                    ps_turb_x = tmpLayer.vX*self.tel.samplingTime
                    ps_turb_y = tmpLayer.vY*self.tel.samplingTime
                    tmpLayer.ratio[0] = ps_turb_x/self.ps_loop
                    tmpLayer.ratio[1] = ps_turb_y/self.ps_loop
                    setattr(self,'layer_'+str(i_layer+1),tmpLayer )
                # self.print_atm()
                
    @property
    def windDirection(self):
         return self._windDirection
    
    @windDirection.setter
    def windDirection(self,val):
        self._windDirection = val

        if self.hasNotBeenInitialized is False:
            if len(val)!= self.nLayer:
                print('Error! Wrong value for the wind-speed! Make sure that you inpute a wind-speed for each layer')
            else:
                print('Updating the wind direction...')
                for i_layer in range(self.nLayer):
                    tmpLayer = getattr(self,'layer_'+str(i_layer+1))
                    # tmpLayer.notDoneOnce = True
                    tmpLayer.direction = val[i_layer]
                    tmpLayer.vY            = tmpLayer.windSpeed*np.cos(np.deg2rad(tmpLayer.direction))                    
                    tmpLayer.vX            = tmpLayer.windSpeed*np.sin(np.deg2rad(tmpLayer.direction))
                    ps_turb_x = tmpLayer.vX*self.tel.samplingTime
                    ps_turb_y = tmpLayer.vY*self.tel.samplingTime
                    tmpLayer.ratio[0] = ps_turb_x/self.ps_loop
                    tmpLayer.ratio[1] = ps_turb_y/self.ps_loop
                    setattr(self,'layer_'+str(i_layer+1),tmpLayer )
                # self.print_atm()


                          
             
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
    def show(self):
        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        print(self.tag+':')
        for a in attributes:
            if not(a[0].startswith('__') and a[0].endswith('__')):
                if not(a[0].startswith('_')):
                    if not np.shape(a[1]):
                        tmp=a[1]
                        try:
                            print('          '+str(a[0])+': '+str(tmp.tag)+' object') 
                        except:
                            print('          '+str(a[0])+': '+str(a[1])) 
                    else:
                        if np.ndim(a[1])>1:
                            print('          '+str(a[0])+': '+str(np.shape(a[1])))  
                            
class Layer:
    pass


        
        
            
            
            
            
    
        