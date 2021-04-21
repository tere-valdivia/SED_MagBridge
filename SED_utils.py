import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt

# Constantes
h = 6.626e-34  # [J seg]
k_b = 1.38e-23  # [J/K]
c = 3.e8  # [m/s]

# funciones


def B_nu(nu, T):
    return ((2 * h * (nu**3)) * (c**(-2))) * ((sp.exp((h * nu) *
                                                      ((k_b * T)**(-1))) - 1)**(-1))
    # planck function [J s^-1 m^-2 sr^-1 Hz^-1]


def B_lam(lam, T):
    return ((2 * h * (c**2)) * (lam**(-5))) * ((sp.exp((h * c) *
                                                       ((lam * k_b * T)**(-1))) - 1)**(-1))
    # planck function [J s^-1 m^-2 sr^-1 m^-1]


def S_nu(nu, cte, beta, T):
    return cte * (nu**beta) * B_nu(nu, T) * 1e26  # [Jy]


def S_lam(lam, cte, beta, T):
    return cte * ((c / lam)**beta) * B_lam(lam, T)


def nuS_nu(nu, cte, beta, T):
    return nu * S_nu(nu, cte, beta, T)  # [Jy Hz]


def chi_square0(entrada):
    '''
    Funcion que toma los flujos reales y la funcion nuS_nu y devuelve el valor
    de chi cuadrado (la suma de las diferencias) luego de utilizar la
    constante ingresada, el factor beta y la temperatura

    Toma entrada = [cte, beta, T]
    '''
    ctes0 = entrada[0]
    betas0 = entrada[1]
    temperaturas0 = entrada[2]
    a = []
    for l in range(len(flujos_fit)):
        a.append((flujos_fit[l] * (frecuencias_fit[l] / 1.e9) - (1.e-9) *
                  nuS_nu(frecuencias_fit[l], ctes0, betas0, temperaturas0))**2 / errores_fit[l]**2)
    # transforma las unidades de frecuencia a GHz, pues estan en Hz
    # parece que lo tiene con normalizacion dependiendo del error que tenga
    return sp.sum(a)


def plot_sed(*arguments):
    """
    Inputs:
    Title (or object)
    Energy (flux * frequency, or y-axis)
    Wavelengths (or x-axis)
    Errors in energy (or y-axis errors)
    Labels of each dot
    Colors of each dot
    Format of each dot (shape)
    Error bars of each dot
    """
    plt.xlabel('$\lambda$ [$\mu$m]')
    plt.ylabel(r'$\nu$S$_{\nu}$ [Jy GHz]')
    source = arguments[0]
    plt.title(str(source))
    # plt.axes.tick_params(length=9)
    plt.yscale('log')
    plt.xscale('log')
    mpl.rcParams['font.size'] = 12.0
    mpl.rcParams['font.family'] = 'serif'
    mpl.rcParams['font.serif'] = 'Times New Roman'
    mpl.rcParams['legend.numpoints'] = 1
    mpl.rcParams['legend.markerscale'] = 1
    mpl.rcParams['legend.handlelength'] = 0
    mpl.rcParams['lines.markersize'] = 8
    mpl.rcParams['xtick.minor.size'] = 0
    mpl.rcParams['ytick.minor.size'] = 0
    mpl.rcParams['xtick.major.size'] = 0
    mpl.rcParams['ytick.major.size'] = 0
    colores = arguments[5]  # RGB
    fluxes = arguments[1]
    eje_x = arguments[2]
    errorsy = arguments[3]
    #errorsx = sp.zeros(len(errorsy))
    leyenda = arguments[4]
    formats = arguments[6]
    errorsabove = arguments[7]
    for i in range(len(fluxes)):
        plt.errorbar(eje_x[i], fluxes[i], yerr=errorsy[i], fmt=formats[i], barsabove=True,
                     color=colores[i], label=leyenda[i], uplims=errorsabove[i], capsize=3)
    plt.legend(frameon=True, loc='lower left', framealpha=0)
    pass
