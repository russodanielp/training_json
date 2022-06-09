#%%
# curveP
import numpy as np

CONCLIST = [1e-6, 2.679636e-06, 4.019455e-06, 6.029182e-06, 9.043773e-06, 1.356566e-05, 2.034849e-05, 3.052273e-05,
            4.57841e-05, 6.867615e-05, 0.0001030142, 0.0001545213, 0.000231782, 0.000347673, 0.0005215095,
            0.0007822643, 0.001173396, 0.001760095, 0.002640142, 0.003960213, 0.005940319, 0.008910479, 0.01336572,
            0.02004858, 0.03007287, 0.0451093, 0.06766395, 0.1014959, 0.1522439, 0.2283658, 0.3425487, 0.5138231,
            0.7707347, 1.156102, 1.734153, 2.601229, 3.901844, 5.852766, 8.77915, 13.16872, 19.75309,
            29.62963, 44.44444, 66.66667, 100]

def curveP(data,
           RNP=100, # max response
           THR=15, # baseline threshold
           MXDV=5, # maximum allowed deviation from monotonicity
            ):
    """ take from alex books chapter curveP https://link.springer.com/content/pdf/10.1007/978-1-0716-2213-1.pdf """


    data = data.sort_values('log(conc)')
    curve_direction = data.iloc[0].resp - data.iloc[-1].resp


    # check for baseline noise
    data.loc[data.resp.abs() < THR, 'resp'] = 0
    # curve direction is > 0
    # means descendning
    # change max response direction
    data.loc[data.resp.abs() > RNP, 'resp'] = -1 * RNP if curve_direction > 0 else RNP

    # correct blips
    if (abs(data.iloc[0].resp) > 0) and (data.iloc[1].resp == 0 and data.iloc[2].resp == 0):
        data.loc[data.index[0], 'resp'] = 0



    # find poiints that violate global directions
    gds = (abs(curve_direction) + MXDV < (data.resp - data.iloc[0].resp).abs())

    if gds.sum() < 2:
        # step 9
        pass

    # get directions
    # if increasing, -1
    # if decreasing, 1
    diffs = data.resp.diff().fillna(0)
    signs = diffs.copy()
    signs[diffs.abs() <= MXDV] = 0
    signs[(diffs.abs() > MXDV) & (diffs > 0)] = -1
    signs[(diffs.abs() > MXDV) & (diffs < 0)] = 1

    # data points not in the correct direction
    correct_direction = curve_direction*signs >= 0

    # assign incorrect data points
    # the previous value
    data_corrected = data.copy()
    data_corrected.loc[~correct_direction, 'resp'] = np.nan
    data_corrected['resp'] = data_corrected['resp'].fillna(method='ffill')
    return data_corrected
