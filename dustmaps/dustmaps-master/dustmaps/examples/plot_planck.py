#!/usr/bin/env python
#
# plot_planck.py
# An example of how to query the Planck Collaboration dust map.
#
# Copyright (C) 2016  Gregory M. Green
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#

from __future__ import print_function, division

import numpy as np
import os.path

try:
    import PIL.Image
except ImportError as error:
    print('This example requires Pillow or PIL.\n'
          'See <http://pillow.readthedocs.io/en/stable/installation.html>.')
    raise error

from astropy.coordinates import SkyCoord
import astropy.units as u

from dustmaps.planck import PlanckQuery


def numpy2pil(a, vmin, vmax):
    a = np.clip((a - vmin) / (vmax - vmin), 0., 1.)
    a = (254.99 * a).astype('u1')
    return PIL.Image.fromarray(a)


def main():
    w,h = (2056,1024)
    l_0 = 0.

    # Create a grid of coordinates
    print('Creating grid of coordinates...')
    l = np.linspace(-180.+l_0, 180.+l_0, 2*w)
    b = np.linspace(-90., 90., 2*h+2)
    b = b[1:-1]
    l,b = np.meshgrid(l, b)

    l += (np.random.random(l.shape) - 0.5) * 360./(2.*w)
    b += (np.random.random(l.shape) - 0.5) * 180./(2.*h)

    coords = SkyCoord(l*u.deg, b*u.deg, frame='galactic')

    planck_components = [
        ('ebv', 0., 1.5),
        ('radiance', 0., 1.5),
        ('tau', 0., 1.5),
        ('temp', 15.*u.K, 25.*u.K),
        ('err_temp', 0.*u.K, 4.*u.K),
        ('beta', 1., 3.),
        ('err_beta', 0., 0.2)]

    for component,vmin,vmax in planck_components:
        # Set up Planck query object
        print('Loading Planck map...')
        planck = PlanckQuery(component=component)

        print('Querying map...')
        res = planck.query(coords)

        # Convert the output array to a PIL image and save
        print('Saving image...')
        img = numpy2pil(res[::-1,::-1], vmin, vmax)
        img = img.resize((w,h), resample=PIL.Image.LANCZOS)
        fname = 'planck_{}.png'.format(component)
        img.save(fname)

    return 0


if __name__ == '__main__':
    main()
