'''
Subclass of Telescope, specifically for Bok
'''

from .telescope import Telescope

class Bok(Telescope):

    def __init__(self, datadir:str, obj_name:str, standard_name:str, debug:bool=False):
        '''
        A class to represent data from the Bok telescope

        Args:
            datadir [str]: The path to the directory with the data organized with the 
                           following subdirectories:
                           1) obj_name (provided as an argument here)
                           2) standard_name (provided as an argument here)
                           3) flats
                           4) zero
                           The obj_name subdirectory should have an obj_name_arc.fits
                           as the arc file for that object
            obj_name [str]: The name of the object subdirectory to reduce
            standard_name [str]: The name of the standard star subdirectory
        '''

        # use the default init method built into Telescope
        super().__init__(datadir, obj_name, standard_name, debug)

        # set some constants just for bok
        # tuning these hyperparameters will change the result!!
        # for extracting the 1D spectrum
        self.MAX_LIGHT_FRAC = 0.1 # maximum light fraction for computing the spectral region
        self.INIT_WIDTH = 1 # initial width for the spectral extraction region
        
        # Parameters for the wavelength fitting
        self.MIN_WAVE = 4700  # for V Filter
        self.MAX_WAVE = 7000  # for V Filter
        self.CALIB_ELEMENTS = ['He', 'Ne', 'Ar']
        
        # pack these into a dictioanry of properties 
        self.PROPS = dict(
            max_frac = self.MAX_LIGHT_FRAC,
            width_init = self.INIT_WIDTH,
            min_wave = self.MIN_WAVE,
            max_wave = self.MAX_WAVE,
            element_list = self.CALIB_ELEMENTS
        )
        
        
    def reduce_data(self, debug=False):
        '''
        Reduce data from the Bok Telescope
        '''

        # the data is read in with the constructor
        # so now we can just start reducing!

        # first just extract the 1d spectrum of the science case
        spec1d = super().extract_1d(self.science_red.data, debug=debug, **self.PROPS)
        stan1d = super().extract_1d(self.stan_red.data, debug=debug, **self.PROPS)
        arc1d = super().extract_1d(self.arc_red.data, arc=True,
                                   background_subtract=False, debug=debug, **self.PROPS)
        
        # then perform the wave solution
        wave = super().wavelength_solver(arc1d, debug=debug, **self.PROPS)

        # finally perform the flux calibration
        flux = super().calibrate_flux(wave, spec1d, stan1d, self.standard_name,
                                      self.standard_dir, debug=debug, **self.PROPS)
        
        
        return wave, flux
    
