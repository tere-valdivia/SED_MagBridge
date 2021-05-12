import scipy as sp
from scipy import optimize as op
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from SED_utils import *
"""
The ALMA point is an upper limit so it is placed by hand
"""
# ff_870mu = 0.117
# ff_1300mu = 0.122
# CO32 = 0.85
# We need to calculate these in the future
# ingreso de datos a partir de un archivo (mJy)
df = pd.read_csv('aphot_MagBridgeF.csv', sep=",",
                 header=0, index_col=0)
# SPT_70mu = df.loc['70']['Flux (mJy)']
HER_100mu = df.loc['images/MagBridgeF_SMC.HERITAGE.PACS100.img_43']['flux'] * 1e3
HER_160mu = df.loc['images/MagBridgeF_SMC.HERITAGE.PACS160.img_43']['flux'] * 1e3
# SPT_160mu = df.loc['160 (Spitzer)']['flux']
HER_250mu = df.loc['images/MagBridgeF_SMC_250um_combined_20121215_img_43']['flux'] * 1e3
HER_350mu = df.loc['images/MagBridgeF_SMC_350um_combined_20121215_img_43']['flux'] * 1e3
HER_500mu = df.loc['images/MagBridgeF_SMC_500um_combined_20121215_img']['flux'] * 1e3
LAB_870mu = df.loc['images/MagBridgeF_MagBri_F_sm18_reprojected_43']['flux'] * 1e3
# LAB_870mu2 = df.loc['870 2']['Flux (mJy)'] - ff_870mu - CO32
# LAB_mean = sp.mean([LAB_870mu, LAB_870mu2])
# ALMA_1300mu = 50.086692   # natural
# ALMA_1300mu = 16.5 - ff_1300mu

# err_SPT_70mu = df.loc['70']['Flux error (mJy)']
err_HER_100mu = df.loc['images/MagBridgeA_100um_43arcsec']['flux_error'] * 1e3
err_HER_160mu = df.loc['images/MagBridgeA_160um_43arcsec']['flux_error'] * 1e3
# err_SPT_160mu = df.loc['160 (Spitzer)']['Flux error (mJy)']
err_HER_250mu = df.loc['images/MagBridgeA_250um_43arcsec']['flux_error'] * 1e3
err_HER_350mu = df.loc['images/MagBridgeA_350um_43arcsec']['flux_error'] * 1e3
err_HER_500mu = df.loc['images/MagBridgeA_500um']['flux_error'] * 1e3
# err_LAB_870mu1 = 0.0575+2.49946864117e-05
err_LAB_870mu = df.loc['images/MagBridgeA_BoA_iter_RM_10_reprojected_43arcsec_cut_79_79']['flux_error'] * 1e3
#err_LAB_870mu2 = df.loc['870 2']['Flux error (mJy)']
#err_LAB_mean = sp.sqrt(err_LAB_870mu**2 + err_LAB_870mu2**2)
# err_ALMA_1300mu = 0.0865937388957644
# err_ALMA_1300mu = 5.5

lambdas = sp.array([100., 160., 250., 350., 500., 870.])
frecuencias = c / (lambdas * 1.e-6)  # [Hz]

flujos = sp.array([HER_100mu, HER_160mu,
                   HER_250mu, HER_350mu, HER_500mu, LAB_870mu]) * 1.e-3  # Jy
errflujos = sp.array([err_HER_100mu, err_HER_160mu, err_HER_250mu, err_HER_350mu,
                      err_HER_500mu, err_LAB_870mu]) * 1.e-3  # el error esta en Jy
erroresy = sp.array([err_HER_100mu, err_HER_160mu, err_HER_250mu, err_HER_350mu,
                     err_HER_500mu, err_LAB_870mu]) * frecuencias * 1.e-9 * 1.e-3  # el error esta en Jy GHz

f = sp.arange(frecuencias[len(frecuencias) - 1] - 300.e9, frecuencias[0] + 100.e11,
              ((frecuencias[0] + 100.e11) - ((len(frecuencias) - 1) - 100.e9)) / 1000.)
lam = (c / f) * 1.e6

# datos para fittear
lambdas_fit = sp.array([100., 160., 250., 350., 500.])
frecuencias_fit = c / (lambdas_fit * 1.e-6)
flujos_fit = sp.array([HER_100mu, HER_160mu, HER_250mu, HER_350mu, HER_500mu]) * 1.e-3
errores_fit = sp.array([err_HER_100mu, err_HER_160mu, err_HER_250mu, err_HER_350mu,
                        err_HER_500mu]) * 1.e-3 * frecuencias_fit
lables_lambda = ['100$\mu$m Herschel',
                 '160$\mu$m Herschel',
                 '250$\mu$m Herschel',
                 '350$\mu$m Herschel',
                 '500$\mu$m Herschel',
                 '870$\mu$m LABOCA']

# La funcion a optimizar es nuS_nu

sp.random.seed(42)


def nuSnufit(xdata, *params):
    cte0, beta0, temp0 = params
    return nuS_nu(xdata, cte0, beta0, temp0)


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
        a.append(((flujos_fit[l]*frecuencias_fit[l])*1.e-9-(1.e-9) *
                  nuS_nu(frecuencias_fit[l], ctes0, betas0, temperaturas0))**2/(errores_fit[l]*1.e-9)**2)
    # transforma las unidades de frecuencia a GHz, pues estan en Hz
    # parece que lo tiene con normalizacion dependiendo del error que tenga
    return sp.sum(a)


x0 = (3.46e-29, 1.4, 22)  # entrada
# optimo0, chi0, niter0, nfcalls0, flag0 = op.fmin(chi_square0, x0, maxiter=None,
#                                                  maxfun=None, full_output=1, disp=0)


optimo0, covmat = op.curve_fit(nuSnufit, frecuencias_fit, flujos_fit *
                               frecuencias_fit, x0, errores_fit, absolute_sigma=True)
C0 = optimo0[0]
B0 = optimo0[1]
T0 = optimo0[2]
chi0 = chi_square0((C0, B0, T0))
errorC0, errorB0, errorT0 = sp.sqrt(sp.diag(covmat))

print("T_0: " + str(T0)+'+-'+str(errorT0))
print("B_0: " + str(B0)+'+-'+str(errorB0))
print("C_0: " + str(C0)+'+-'+str(errorC0))
print("Chi0: " + str(chi0))
print("A 870 micrones esperamos (mJy): " +
      str(round(S_nu(frecuencias[5], C0, B0, T0)*1.e3, 2)))
print("A 1.3mm esperamos (mJy): "+str(round(S_nu(frecuencias[6], C0, B0, T0)*1.e3, 2)))
# Guardar los Valores
df = pd.DataFrame()
df['Params'] = ['C', 'Beta', 'T', 'chi_2']
df['Values'] = [C0, B0, T0, chi0]
df['Errors'] = [errorC0, errorB0, errorT0, sp.nan]
df.to_csv('parametros_fit_curvefit_F.csv')

# Plots
fig = plt.figure(figsize=(12, 5))

colores = [(1, 0, 0), (1, 0, 0), (0, 0, 1), (1, 0, 0),
           (1, 0, 0), (0. / 255, 0. / 249, 255. / 255)]  # RGB
fmt = ['.', '.', '.', '.', '.', '.', '_']
errorabove = sp.array([0, 0, 0, 0, 0, 0], dtype=bool)
# Ploteo normal
ax = fig.add_subplot(1, 2, 1)
# los convierte a JyGHz
plt.plot(lam, nuS_nu(f, C0, B0, T0) * 1.e-9, linewidth=1, color='k',
         label=r'S$_{\nu} \propto \nu^{\beta}B_{\nu}(T)$')
plt.plot(lambdas, sp.zeros(len(lambdas)), 'b-')
plt.title('MagBridgeA SED')

# plot_sed('MagBridgeA SED', flujos * (frecuencias / 1.e9),
#          lambdas, erroresy, lables_lambda, colores, fmt, errorabove)
plt.xlabel('$\lambda$ [$\mu$m]')
plt.ylabel(r'$\nu$S$_{\nu}$ [Jy GHz]')
plt.yscale('log')
plt.xscale('log')
# plt.errorbar(lambdas[[2]], flujos[[2]] * (frecuencias[[2]] / 1.e9), yerr=erroresy[[2]],
#              fmt='.', barsabove=True, color=(0, 0, 1), label='Spitzer', capsize=3)
plt.errorbar(lambdas[:5], flujos[:5] * (frecuencias[:5] / 1.e9),
             yerr=erroresy[:5], fmt='.', barsabove=True,
             color=(1, 0, 0), label='Herschel', capsize=3)
plt.errorbar(lambdas[5], flujos[5] * (frecuencias[5] / 1.e9), yerr=erroresy[5],
             fmt=fmt[5], barsabove=True, color=colores[5], label='LABOCA', capsize=3)
# plt.errorbar(lambdas[8], flujos[8] * (frecuencias[8] / 1.e9), yerr=erroresy[8],
#              fmt=fmt[8], barsabove=True, color=colores[8], label='LABOCA n2', capsize=3)
# Solo si se quiere el promedio
# plt.errorbar(lambdas[8], LAB_mean*1.e-3 * (frecuencias[7] / 1.e9), yerr=err_LAB_mean*frecuencias[8]*1.e-12,
#              fmt=fmt[8], barsabove=True, color='gray', label='LABOCA mean', capsize=3)

# plt.errorbar(lambdas[6], flujos[6] * (frecuencias[6] / 1.e9), yerr=erroresy[6],
#              fmt=fmt[6], barsabove=True, color=colores[6], label='ALMA', uplims=True, capsize=3)
plt.figtext(0.4, 0.6, r'$\chi$$^2$$_{red}$ = ' +
            str(round(chi0 / (len(flujos_fit) - 3), 2)) +
            '\nC = ' + str(round(C0, 33)) +
            '\n' + r'$\beta$ = ' + str(round(B0, 1)) + '\nT = '
            + str(round(T0, 1)) + ' K')
plt.legend(frameon=True, loc='lower left', framealpha=0)

plt.ylim(1.e-0, 1.e4)
plt.xlim(6.e1, 1.2e4)

# Excesos
valores = flujos * frecuencias
# print(valores / 1.e9)
modelo = nuS_nu(frecuencias, C0, B0, T0)

excesos = valores / modelo
print(excesos)
ax2 = fig.add_subplot(1, 2, 2)

plt.xlabel('$\lambda$ [$\mu$m]')
plt.ylabel('Residuals (data/model)')
plt.title('MagBridgeA SED Residual')

errorsx = [0, 0, 0, 0, 0, 0]
errorsy = erroresy / (nuS_nu(frecuencias, C0, B0, T0) * 1.e-9)
print(errorsy)
errorabove = sp.array([0, 0, 0, 0, 0, 0], dtype=bool)

# mean_residual = LAB_mean*1.e-3 * frecuencias[8] / modelo[8]
# errory_mean = err_LAB_mean*frecuencias[8]*1.e-12 / (nuS_nu(frecuencias[8], C0, B0, T0) * 1.e-9)
# print('El residual del promedio de LABOCA es ' +
#       str(object=round(mean_residual, 2))+'+-'+str(errory_mean))
residuals = (flujos * frecuencias * 1.e-9) / \
    (nuS_nu(frecuencias, C0, B0, T0) * 1.e-9)

print(residuals)
# Para guardar los residuales en un archivo
# residualTable = pd.DataFrame(data=sp.transpose(
#     [sp.around(residuals, 2), sp.around(errorsy, 2)]), index=lambdas, columns=['E(870)', 'error'])
# residualTable.to_csv('excesos_MagBridgeF_refcorrected.csv')

i = 0

for i in range(len(residuals)):
    plt.errorbar(lambdas[i], residuals[i], errorsy[i], errorsx[i], fmt[i],
                 barsabove=True, markersize=10, capsize=3, mec='k', mew=0.8,
                 elinewidth=2, linewidth=2, color=colores[i], ecolor=colores[i],
                 label=lables_lambda[i], uplims=errorabove[i])
# plt.errorbar(lambdas[8], mean_residual, yerr=errory_mean, fmt='.', barsabove=True, markersize=10, capsize=3, mec='k', mew=0.8,
#              elinewidth=2, linewidth=2, color='gray', ecolor='gray')
line = sp.arange(0, 1500, 10)
plt.xlim(50, 1500)
plt.plot(line, sp.ones(len(line)), 'k-')
# plt.savefig('MagBridgeA_SED_nov2019_refcorrected.eps',
#             format='eps', transparent=True)
# plt.savefig('MagBridgeA_SED_jul2019_withLABmean.eps',
#             format='eps', transparent=True)
plt.show()


# def masscont(fluxes, nu, epsilon, Td, distance, mu=1.36):
#     flux, fluxerror = fluxes
#     eps, epserr = epsilon
#     dist = distance / 3.24078e-19
#     Bd = B_nu(nu, Td) * 1.e26  # Jy sr-1
#     mh = 1.6e-27
#     sunmass = 1.9e30
#     mass = flux * mu * mh * (dist)**2 / (Bd * eps) / sunmass
#     masserror = sp.sqrt((epserr/eps)**2+(fluxerror/flux)**2) * mass
#     return sp.array([mass, masserror])
#
#
# masscont(sp.array([LAB_870mu, err_LAB_870mu])*1.e-3,
#          frecuencias[6], sp.array([3.94, 0.05])*1.e-27, T0, 60000)
# modelo = nuS_nu(frecuencias, C0, B0, T0)
# valorespredichos = S_nu(frecuencias, C0, B0, T0)
# e_valorespredichos = sp.sqrt((errorT0/T0)**2+(errorB0/B0)**2) * valorespredichos
#
# kappa160 = sp.array([9.6, 2.9])
# epsilon160 = kappa160 * 0.0015 * 1.36 * 1.6e-24
# masses = masscont((valorespredichos[1], e_valorespredichos[1]), frecuencias[1],
#                   epsilon160, T0, 60000)
