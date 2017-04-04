"""
Script to generate sky maps
with a powerlaw power spectrum
"""

import numpy
from numpy import pi
import healpy
import matplotlib.pyplot as plt

lognside = 9
nside = 2**lognside
npix = healpy.pixelfunc.nside2npix(nside)
corrsize = 3 * nside
l = numpy.arange(corrsize)
# convert to angular distance in degrees
deg = 2*pi / (1+l) * 180 / pi
pixsize = healpy.nside2pixarea(nside, degrees=True)**0.5
print 'pixel size:', pixsize
print 'degrees:'
print deg
print 'smallest: %.4f arcmin' % (deg.min() * 60)

if False:
	# get power spectrum from real map
	values = numpy.random.normal(size=npix)
	plt.figure()
	healpy.visufunc.mollview(values)
	plt.savefig('genucorr.pdf', bbox_inches='tight')
	plt.savefig('genucorr.png', bbox_inches='tight')
	plt.close()
	r = healpy.sphtfunc.anafast(values)
	assert len(r) == corrsize
	plt.plot(deg, r)

# given a power spectrum, generate a new map
# use power spectrum slope of Markowitz+14
# it is 1.5, but measured in days. Need to assume
# something about the orbits to convert to angles
# E.g. keplerian circular orbits give only a 
# scale conversion factor
r = 1e-3 * (deg/(2*180))**0.8 + 1e-10
r[0] = 1e-10
assert numpy.isfinite(r).all(), r
plt.plot(deg, r)
plt.xscale('log')
plt.yscale('log')
plt.savefig('gen_corrfuncs.pdf')
plt.savefig('gen_corrfuncs.png')
plt.close()
print 'generating...'
values = healpy.sphtfunc.synfast(r, nside)
print 'generating done'
assert numpy.isfinite(values).all(), values
plt.figure()
healpy.visufunc.mollview((values < 0)*1, cmap='gray_r')
plt.title('Obscured fraction: %.2f%%' % ((values < 0).mean()*100))
plt.savefig('gencorr.pdf', bbox_inches='tight')
plt.savefig('gencorr.png', bbox_inches='tight')
plt.close()



