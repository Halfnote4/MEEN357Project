import subfunctions
import numpy

Crr_array = numpy.linspace(0.01,0.4,25)
slope_array_deg = numpy.linspace(-10,35,25)



CRR, SLOPE = numpy.meshgrid(Crr_array, slope_array_deg)

VMAX = numpy.zeros(numpy.shape(CRR), dtype = float)


N = numpy.shape(CRR)[0]
for i in range(N):
    for j in range(N):
        Crr_sample = float(CRR[i,j])
        slope_sample = float(SLOPES[i,j])
        VMAX[i,j] = ... # here you put code to find the max speed at Crr_sample and
 # slope_sample
