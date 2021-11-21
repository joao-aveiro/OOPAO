import inspect

import numpy as np

from photometry import Photometry


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CLASS INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        _ src.nPhoton   : number of photons per m2 per s. if this property is changed after the initialization,
                            the magnitude is automatically updated to the right value.
        _ src.fluxMap   : 2D map of the number of photons per pixel per frame (depends on the loop frequency
                            defined by tel.samplingTime)
                 
        ************************** EXAMPLE **************************

        Create a source object in H band with a magnitude 8 and combine it to the telescope
        src = Source(opticalBand = 'H', magnitude = 8) 
        src*tel

        
        """
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
        tmp = self.photometry(opt_band)                         # get the photometry properties
        self.opt_band = opt_band                                # optical band
        self.wavelength = tmp["wl"]                             # wavelength in m
        self.bandwidth = tmp["bw"]                              # optical bandwidth
        self.zero_point = tmp["zp"]                             # zero point
        self.magnitude = magnitude                              # magnitude
        self.phase = []                                         # phase of the source
        self.phase_no_pupil = []                                # phase of the source
        self.flux_map = []                                      # 2D flux map of the source
        self.n_photon = self.zero_point*10**(-0.4*magnitude)    # number of photon per m2 per s
        self.tag = 'source'                                     # tag of the object

        print(f"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SOURCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
              f"Wavelength \t{round(self.wavelength*1e6, 3)} \t [microns]\n"
              f"Optical Band \t{opt_band}\n"
              f"Magnitude \t{self.magnitude}\n"
              f"Flux \t\t{np.round(self.n_photon)}\t [photons/m2/s]\n"
              '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SOURCE INTERACTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    def __mul__(self, tel_object):
        # update the phase of the source
        self.phase = tel_object.OPD * 2 * np.pi / self.wavelength
        self.phase_no_pupil = tel_object.OPD_no_pupil * 2 * np.pi / self.wavelength

        # compute the variance in the pupil
        self.var = np.var(self.phase[np.where(tel_object.pupil == 1)])
        # assign the source object to the telescope object
        tel_object.src = self

        self.fluxMap = tel_object.pupilReflectivity * self.n_photon*tel_object.samplingTime\
                       * (tel_object.D/tel_object.resolution)**2
        
        return tel_object

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SOURCE PHOTOMETRY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    @staticmethod
    def photometry(arg):
        # photometry object [wavelength, bandwidth, zeroPoint]
        phot = Photometry()
        if isinstance(arg, str):
            if hasattr(phot, arg):
                return getattr(phot, arg)
            else:
                print("Error: Wrong name for the photometry object.")
                return -1
        else:
            print("Error: The photometry object takes a scalar as an input.")
            return -1

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SOURCE PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    @property
    def n_photon(self):
        return self._n_photon
    
    @n_photon.setter
    def n_photon(self, val):
        self._n_photon = val
        self.magnitude = -2.5*np.log10(val/self.zero_point)
        print("NGS flux updated!"
              "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SOURCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
              f"Wavelength \t{round(self.wavelength*1e6, 3)}\t [microns]"
              f"Optical Band \t{self.opt_band}"
              f"Magnitude \t{self.magnitude}"
              f"Flux \t\t{np.round(self.n_photon)}\t [photons/m2/s]')"
              "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

#    @property
#    def magnitude(self):
#        return self._magnitude
#    
#    @magnitude.setter
#    def magnitude(self,val):
#        self._magnitude  = val
#        self.nPhoton     = self.zeroPoint*10**(-0.4*val)
#            
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    def show(self):
        attributes = inspect.getmembers(self, lambda a: not(inspect.isroutine(a)))
        print(self.tag+':')
        for a in attributes:
            if not(a[0].startswith('__') and a[0].endswith('__')):
                if not(a[0].startswith('_')):
                    if not np.shape(a[1]):
                        tmp = a[1]
                        try:
                            print('          '+str(a[0])+': '+str(tmp.tag)+' object') 
                        except:
                            print('          '+str(a[0])+': '+str(a[1])) 
                    else:
                        if np.ndim(a[1]) > 1:
                            print('          '+str(a[0])+': '+str(np.shape(a[1])))   
