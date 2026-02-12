'''
Miscellaneous utilities for LIM code
'''
import numpy as np
from astropy.units.quantity import Quantity
from astropy.cosmology.core import FlatLambdaCDM
from scipy.interpolate import interp1d
import inspect
import astropy.units as u

import luminosity_functions as lf
import mass_luminosity as ml

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors

class cached_property(object):
    """
    From github.com/Django, who wrote a much better version of this than
    the one I had previously.
    
    Decorator that converts a self.func with a single self argument into a
    property cached on the instance.
    """
    def __init__(self, func):
        self.func = func

    def __get__(self, instance, type=None):
        if instance is None:
            return self
        
        # ADDED THIS CODE TO LIST PROPERTY FOR UPDATING
        instance._update_list.append(self.func.__name__)
        
        res = instance.__dict__[self.func.__name__] = self.func(instance)
        return res
        
def get_default_params(func):
    '''
    Gets the default parameters of a function as input to check_params. Output
    is a dictionary of parameter names and values. If the function has a
    "self" argument it is removed from the dictionary.
    '''
    
    args = inspect.getargspec(func)
    
    param_names = args.args
    if 'self' in param_names:
        param_names.remove('self')
    
    default_values = args.defaults
    
    default_params = dict(zip(param_names,default_values))

    return default_params
    
        
def check_params(input_params, default_params):
    '''
    Check input parameter values to ensure that they have the same type and
    unit as the required inputs
    '''
    
    for key in input_params.keys():
        # Check if input is a valid parameter
        if key not in default_params.keys():
            raise AttributeError(key+" is not a valid parameter")
        
        input_value = input_params[key]
        default_value = default_params[key]
        
        # Check if input has the correct type
        if type(input_value)!=type(default_value):
            # Some inputs can have multiple types
            if key=='cosmo_model':
                if type(input_value)==FlatLambdaCDM:
                    pass
                else:
                    raise(TypeError(
                      "Parameter cosmo_model must be a str or FlatLambdaCDM"))
            elif key=='scatter_seed':
                if type(input_value)==int or type(input_value)==float:
                    pass
            
            elif key=='Nfeeds':
                if type(input_value)==int or type(input_value)==float:
                    pass
                
            elif type(default_value)==Quantity:
                raise TypeError("Parameter "+key+
                        " must be an astropy quantity")
            else:
                raise TypeError("Parameter "+key+" must be a "+
                                    str(type(default_value)))
            
        # If input is a quantity, check if it has the correct dimension
        elif (type(default_value)==Quantity and not
                 input_value.unit.is_equivalent(default_value.unit)):
            
            # Tmin/Tmax may be in either uK or Jy/sr depending on do_Jysr     
            if key=='Tmin' or key=='Tmax':
                if (input_params['do_Jysr'] and 
                   input_value.unit.is_equivalent(u.Jy/u.sr)):
                    pass
                else:
                    raise TypeError("Parameter "+key+
                                " must have units equivalent to "
                                +str(default_value.unit))
                                
        # Special requirements for certain parameters
        elif (key=='model_type' and not 
                (input_value=='ML' or input_value=='LF')):
            # model_type can only be ML or LF
            raise ValueError("model_type must be either 'ML' or 'LF'")
            
            
def check_model(model_type,model_name):
    '''
    Check if model given by model_name exists in the given model_type
    '''
    if model_type=='ML' and not hasattr(ml,model_name):
        if hasattr(lf,model_name):
            raise ValueError(model_name+" not found in mass_luminosity.py."+
                    " Set model_type='LF' to use "+model_name)
        else:
            raise ValueError(model_name+
                    " not found in mass_luminosity.py")
    elif model_type=='LF' and not hasattr(lf,model_name):
        if hasattr(ml,model_name):
            raise ValueError(model_name+
                    " not found in luminosity_functions.py."+
                    " Set model_type='ML' to use "+model_name)
        else:
            raise ValueError(model_name+
                    " not found in luminosity_functions.py")
            
                                
def ulogspace(xmin,xmax,nx):
    '''
    Computes logarithmically-spaced numpy array between xmin and xmax with nx
    points.  This function calls the usual np.loglog but takes the linear
    values of xmin and xmax (where np.loglog requires the log of those limits)
    and allows the limits to have astropy units.  The output array will be a
    quantity with the same units as xmin and xmax
    '''
    xmax = xmax.to(xmin.unit)
    return np.logspace(np.log10(xmin.value),np.log10(xmax.value),nx)*xmin.unit
    
def ulinspace(xmin,xmax,nx):
    '''
    Computes linearly-spaced numpy array while preserving astropy units
    '''
    xmax = xmax.to(xmin.unit)
    return np.linspace(xmin.value,xmax.value,nx)*xmin.unit
    
def uinterp1d(x,y,**kwargs):
    '''
    A version of the scipy interpolator that keeps astropy units on the output
    value.  Should accept the same inputs as scipy.interp1d, but can accept
    dimensionful x and y and outputs an array with the same units as y.
    '''
        
    if hasattr(y,'unit'):
        yunit = y.unit
        y = y.value
    else:
        yunit = 1.0
        
    if hasattr(x,'unit'):
        xunit = x.unit
        x = x.value
        return lambda xi: interp1d(x,y,**kwargs)(xi.to(xunit).value)*yunit
        
    else:
        return lambda xi: interp1d(x,y,**kwargs)(xi)*yunit
    
def dictCombine_2(x,y):
    '''
    Return a single dict containing all the elements of x and y
    '''
    # Check if x and y have common elements
    for key in y:
        if key in x:
            print("Warning: Overlapping keys. Value in y will be used")
            
    return dict(x,**y)
    
def dictCombine_3(x,y,z):
    '''
    Returns a single dict containing all the elements of x, y, and z
    '''
    for key in z:
        if (key in y) or (key in x):
            print("Warning: Overlapping keys")
            
    for key in y:
        if key in x:
            print("Warning: Overlapping keys")
    
    return dict(x,**dict(y,**z))
    
def LtoLprime(L,nu):
    '''
    Convert luminosity in solar luminosities to K km/s pc^2 units
    '''
    L0 = L.to(u.Lsun).value
    nu0 = nu.to(u.GHz).value
    return L0/4.9e-5*(nu0/115)**-3 *u.K*u.km*u.s**-1*u.pc**2
    
    
    
class MidpointNormalize(colors.Normalize):
	"""
	Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)

	e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
	"""
	def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
		self.midpoint = midpoint
		colors.Normalize.__init__(self, vmin, vmax, clip)

	def __call__(self, value, clip=None):
		# I'm ignoring masked values and all kinds of edge cases to make a
		# simple example...
		x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
		return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))
    
def plotCov(m,vmin=-1,vmax=1):
    '''
    Plot covariance matrix with bwr colormap with uncorrelated points
    colored white
    '''
    m[np.isnan(m)] = 0.
    plt.imshow(m,cmap=matplotlib.cm.bwr,norm=MidpointNormalize(midpoint=0,vmin=vmin,vmax=vmax),extent=(-2,1,1,-2))
    
def cov_to_corr(covmat):
        '''
        Converts covariance matrix to correlation matrix (with unity on the
        diagonal)
        '''
        s = np.sqrt(covmat.diagonal())
        s1,s2 = np.meshgrid(s,s)
        return covmat/(s1*s2)
    
def conc22(m1,m2,m3,m4):
    '''
    Concatenate four square matrices into one big one
    '''
    return np.concatenate((np.concatenate((m1,m2)),
                np.concatenate((m3,m4))),axis=1)
                
################
# Rebin arrays #
################

def rebin_map(pix,rebin,operation='mean',dim_const=0):
    '''
    Rebin 3d intensity map in the plane of the sky by averaging together rebin 
    x rebin groups of pixels
    '''
    c,a,b = pix.shape
    pix = pix[:,0:a-a%rebin,0:b-b%rebin]

    c1,a1,b1 = pix.shape

    if dim_const==0:
        return bin_ndarray(pix,(c1,a1/rebin,b1/rebin),operation=operation)
    elif dim_const==1:
        return bin_ndarray(pix,(c1/rebin,a1,b1/rebin),operation=operation)
    else:
        return bin_ndarray(pix,(c1/rebin,a1/rebin,b1),operation=operation)
        

def bin_ndarray(ndarray, new_shape, operation='mean'):
    """
    Bins an ndarray in all axes based on the target shape, by summing or
        averaging.

    Number of output dimensions must match number of input dimensions and 
        new axes must divide old ones.

    Example
    -------
    >>> m = np.arange(0,100,1).reshape((10,10))
    >>> n = bin_ndarray(m, new_shape=(5,5), operation='sum')
    >>> print(n)

    [[ 22  30  38  46  54]
     [102 110 118 126 134]
     [182 190 198 206 214]
     [262 270 278 286 294]
     [342 350 358 366 374]]

    """
    operation = operation.lower()
    if not operation in ['sum', 'mean']:
        raise ValueError("Operation not supported.")
    if ndarray.ndim != len(new_shape):
        raise ValueError("Shape mismatch: {} -> {}".format(ndarray.shape,
                                                           new_shape))
    compression_pairs = [(d, c//d) for d,c in zip(new_shape,
                                                  ndarray.shape)]
    flattened = [l for p in compression_pairs for l in p]
    ndarray = ndarray.reshape(flattened)
    for i in range(len(new_shape)):
        op = getattr(ndarray, operation)
        ndarray = op(-1*(i+1))
    return ndarray             

#################
# Wiener filter #
#################
class WF(object):
    '''
    Class for performing basic Wiener filtering
    '''
    def __init__(self,d,N,S=None):
        self.d = d
        self._S = S
        self.N = N
        
        self._update_list = []
        
    #@cached_property
    #def d(self):
    #    '''
    #    Data vector, make sure it's a 1xN array
    #    '''
    #    if self._d.ndim==1:
    #        return self._d.reshape(1,self._d.size)
    #    elif self._d.ndim==2:
    #        return self._d
        
    @cached_property
    def S(self):
        if self._S is None:
            return np.diag(self.d**2)
        else:
            return self._S
            
    @cached_property
    def SN_inv(self):
        '''
        Inverse of S+N
        '''
        return np.linalg.inv(self.S+self.N)
    
    @cached_property
    def inv_check(self):
        '''
        Check if S+N inversion is stable
        '''
        return np.matmul(self.SN_inv,self.S+self.N)
        
    @cached_property
    def F(self):
        '''
        Wiener filter matrix
        '''
        return np.matmul(self.S,self.SN_inv)
        
    @cached_property
    def s(self):
        '''
        Filtered signal
        '''
        dv = np.reshape(self.d,(1,self.d.size))
        return np.matmul(self.F,dv.T)
        
    @cached_property
    def C_s(self):
        '''
        Covariance of filtered signal
        '''
        return np.matmul(self.F,self.N)
        
    @cached_property
    def sig_s(self):
        '''
        Error on filtered signal
        '''
        return np.sqrt(np.diag(self.C_s))
