#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from stvid.stio import fourframe
import ephem
from astropy.time import Time
import datetime
from astropy.coordinates import SkyCoord, FK5
from astropy.wcs import WCS
import astropy.units as u
import glob

class two_line_element:
    """TLE class"""

    def __init__(self, tle0, tle1, tle2):
        """Define a tle"""

        self.tle0 = tle0 #.decode("utf-8")
        self.tle1 = tle1 #.decode("utf-8")
        self.tle2 = tle2 #.decode("utf-8")
        if self.tle0[:2]=="0 ":
            self.name = self.tle0[2:].strip()
        else:
            self.name = self.tle0.strip()
        self.id = self.tle1.split(" ")[1][:5]

# Select those TLEs inside the fov
def select_tles(tles, observer, t, pcen, radius):
    # Loop over all tles
    ra = []
    dec = []
    for tle in tles:
        # Set satellite
        satellite = ephem.readtle(str(tle.tle0),
                                  str(tle.tle1),
                                  str(tle.tle2))
        
        # Set observer date
        observer.date = ephem.date(t.datetime)

        # Compute satellite observer
        satellite.compute(observer)
        ra.append(float(satellite.ra))
        dec.append(float(satellite.dec))

    # Convert to skycoords
    pos = SkyCoord(ra=ra, dec=dec, unit="rad", frame="fk5", equinox=tstart).transform_to(FK5(equinox="J2000"))

    # Compute separations
    r = pos.separation(pcen)
    c = r<radius

    return tles[c]
    
def plot_objects(tles, observer, t):
    # Loop over all tles
    for tle in tles:
        # Set satellite
        satellite = ephem.readtle(str(tle.tle0),
                                  str(tle.tle1),
                                  str(tle.tle2))

        # Loop over times
        ra = []
        dec = []
        for t0 in t:
            # Set observer date
            observer.date = ephem.date(t0.datetime)

            # Compute satellite observer
            satellite.compute(observer)
            ra.append(float(satellite.ra))
            dec.append(float(satellite.dec))

        # Convert to skycoords
        pos = SkyCoord(ra=ra, dec=dec, unit="rad", frame="fk5", equinox=tstart).transform_to(FK5(equinox="J2000"))

        ax.plot(pos[0].ra, pos[0].dec, transform=ax.get_transform("world"), marker=".", color="k")
        ax.plot(pos.ra, pos.dec, transform=ax.get_transform("world"), color="k")
#        ax.plot(pos[0].ra, pos[0].dec, "text", transform=ax.get_transform("world"))
        ax.text(pos[0].ra.degree, pos[0].dec.degree, " %s"%tle.id,
                transform=ax.get_transform("world"),
                color="black", fontsize=10)
        
if __name__ == "__main__":

    # Set observer
    observer = ephem.Observer()
    observer.lon = "6.3785"
    observer.lat = "52.8344"
    observer.elevation = 10

    # Read TLEs
    with open("/data2/tle/catalog.tle", "r") as fp:
        lines = fp.readlines()
    alltles = np.array([two_line_element(lines[i], lines[i+1], lines[i+2]) for i in range(0, len(lines), 3)])


    # Loop over files
    fnames = sorted(glob.glob("/home/bassa/satobs/20190503_205547_4/processed/2019-05-04T02:36:*.fits"))
    for fname in fnames:
        # Read fourframe
        ff = fourframe(fname)
    
        # Dates
        tstart = Time(ff.mjd, format="mjd")
        tmid = tstart+ff.texp*u.s
        t = tstart + np.linspace(0, ff.texp, 10)*u.s
    
        # Center
        pcen = SkyCoord(ra=ff.crval[0], dec=ff.crval[1], unit="deg", frame="icrs")

        # Down select TLEs
        tles = select_tles(alltles, observer, tmid, pcen, 10*u.deg)
    
        # Get wcs
        ff.header["NAXIS"] = 2
        w = WCS(ff.header)
        
        # Image stats
        v = ff.zmax
        vmean, vstd = np.mean(v), np.std(v)
        vmin, vmax = vmean-2.0*vstd, vmean+6.0*vstd

        # Start figure
        fig = plt.figure(figsize=(8, 6.4))
        ax = fig.add_subplot(111, projection=w)

        # Plot images
        ax.imshow(v, origin="lower", aspect=1, interpolation="None", cmap="magma",
                  vmin=vmin, vmax=vmax)
        ax.set_xlabel("Right Ascension")
        ax.set_ylabel("Declination")
        ax.grid(alpha=0.5)
        ax.set_xlim(0, ff.nx)
        ax.set_ylim(0, ff.ny)
    
        plot_objects(tles, observer, t)

        plt.savefig("satid.png")
