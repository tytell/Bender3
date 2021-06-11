import numpy as np
from scipy import interpolate

def make_motor_stepper_pulses(t, pos, vel, outsampfreq,
                              signconvention='Left is positive',
                              scale=6.0,
                              stepsperrev=6400.0):

    tout = np.arange(t[0], t[-1], 1.0/outsampfreq)

    poshi = interpolate.interp1d(t, pos, kind='linear', assume_sorted=True, bounds_error=False,
                                 fill_value=0.0)(tout)
    velhi = interpolate.interp1d(t, vel, kind='linear', assume_sorted=True, bounds_error=False,
                                 fill_value=0.0)(tout)

    if signconvention == 'Left is negative' or \
            signconvention == 'Right is positive':
        poshi = -poshi
        velhi = -velhi
    elif signconvention == 'Left is positive' or \
            signconvention == 'Right is negative':
        pass
    else:
        raise ValueError('Unrecognized motor sign convention {}'.format(signconvention))

    poshi *= scale
    velhi *= scale

    if outsampfreq == 0 or stepsperrev == 0:
        raise ValueError('Problems with parameters')

    stepsize = 360.0 / stepsperrev
    maxspeed = stepsize * outsampfreq / 2

    if np.any(np.abs(vel) > maxspeed):
        raise ValueError('Motion is too fast!')

    stepnum = np.floor(poshi / stepsize)
    dstep = np.diff(stepnum)
    motorstep = np.concatenate((np.array([0], dtype='uint8'), (dstep != 0).astype('uint8')))
    motordirection = (velhi <= 0).astype('uint8')

    motorenable = np.ones_like(motordirection, dtype='uint8')
    motorenable[-5:] = 0

    dig = np.packbits(np.column_stack((np.zeros((len(motorstep), 5), dtype=np.uint8),
                                        motorenable,
                                        motorstep,
                                        motordirection)))

    return tout, dig, motorstep, motordirection
