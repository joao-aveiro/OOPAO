import inspect

import numpy as np

from photometry import Phot


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CLASS INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

class Source:    
    def __init__(self, opt_band, magnitude):
        """
        ************************** REQUIRED PARAMETERS **************************
        
        A Source object is characterised by two parameter:
        _ optBand               : the optical band of the source (see the method photometry)
        _ magnitude             : The magnitude of the star
                            
        ************************** COUPLING A SOURCE OBJECT **************************
        
        Once generated, a Source object "src" can be coupled to a Telescope "tel" that contains the OPD.
        _ This is achieved using the * operator     : src*tel
        _ It can be accessed using                  : tel.src       

    
        ************************** MAIN PROPERTIES **************************
        
        The main properties of a Telescope object are listed here: 
        _ src.phase     : 2D map of the phase scaled to the src wavelength corresponding to tel.OPD
        _ src.nPhoton   : number of photons per m2 per s. if this property is changed after the initialization, the magnitude is automatically updated to the right value. 
        _ src.fluxMap   : 2D map of the number of photons per pixel per frame (depends on the loop frequency defined by tel.samplingTime)  
                 
        ************************** EXEMPLE **************************

        Create a source object in H band with a magnitude 8 and combine it to the telescope
        src = Source(opticalBand = 'H', magnitude = 8) 
        src*tel

        
        """
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  
        tmp             = self.photometry(opt_band)              # get the photometry properties
        self.optBand    = opt_band
        self.wavelength = tmp[0]                                # wavelength in m
        self.bandwidth  = tmp[1]                                # optical bandwidth
        self.zeroPoint  = tmp[2]/368                            # zero point
        self.magnitude  = magnitude                             # magnitude
        self.phase      = []                                    # phase of the source 
        self.phase_no_pupil      = []                                    # phase of the source 
        self.fluxMap    = []                                    # 2D flux map of the source
        self.nPhoton    = self.zeroPoint*10**(-0.4*magnitude)   # number of photon per m2 per s
        self.tag        = 'source'                              # tag of the object
        
        
        print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SOURCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        print('Wavelength \t'+str(round(self.wavelength*1e6,3)) + ' \t [microns]') 
        print('Optical Band \t'+opt_band)
        print('Magnitude \t' + str(self.magnitude))
        print('Flux \t\t'+ str(np.round(self.nPhoton)) + str('\t [photons/m2/s]'))
        print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SOURCE INTERACTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        
    def __mul__(self,telObject):
        # update the phase of the source
        self.phase      = telObject.OPD*2*np.pi/self.wavelength
        self.phase_no_pupil      = telObject.OPD_no_pupil*2*np.pi/self.wavelength

        # compute the variance in the pupil
        self.var        = np.var(self.phase[np.where(telObject.pupil==1)])
        # assign the source object to the telescope object
        telObject.src   = self

        self.fluxMap    = telObject.pupilReflectivity*self.nPhoton*telObject.samplingTime*(telObject.D/telObject.resolution)**2
        
        return telObject
     
    
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SOURCE PHOTOMETRY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    def photometry(self,arg):
        # photometry object [wavelength, bandwidth, zeroPoint]
        phot = Phot()
        if isinstance(arg,str):
            if hasattr(phot,arg):
                return getattr(phot,arg)
            else:
                print('Error: Wrong name for the photometry object')
                return -1
        else:
            print('Error: The photometry object takes a scalar as an input')
            return -1   
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SOURCE PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    @property
    def nPhoton(self):
        return self._nPhoton
    
    @nPhoton.setter
    def nPhoton(self,val):
        self._nPhoton  = val
        self.magnitude = -2.5*np.log10(val/self.zeroPoint)
        print('NGS flux updated!')
        print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SOURCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        print('Wavelength \t'+str(round(self.wavelength*1e6,3)) + ' \t [microns]') 
        print('Optical Band \t'+self.optBand) 
        print('Magnitude \t' + str(self.magnitude))
        print('Flux \t\t'+ str(np.round(self.nPhoton)) + str('\t [photons/m2/s]'))
        print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        
#    @property
#    def magnitude(self):
#        return self._magnitude
#    
#    @magnitude.setter
#    def magnitude(self,val):
#        self._magnitude  = val
#        self.nPhoton     = self.zeroPoint*10**(-0.4*val)
#            
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
            
