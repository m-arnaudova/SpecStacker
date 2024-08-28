# SpecStacker: The Spectral Stacking Code
This is a new rest-frame spectral stacking code for galaxies and AGN spectra. For a quick start-up tutorial, see [example](https://nbviewer.org/github/m-arnaudova/SpecStacker/blob/main/stack_example.ipynb) notebook.

The code takes arrays of observed-frame wavelength, flux density, flux uncertainties and the corresponding redshifts as input parameters, shifts all spectra to the rest-frame, and resamples them onto a common wavelegth grid based on the redshift distribution and an input pixel size (e.g. 1Å). The normalisation of each spectrum is performed by taking the median, computed at the reddest possible end of the area where all spectra populate the common wavelength grid, with prominent emission lines masked out. The composite spectrum is then created by taking the median of all normalized flux density values that fall within a given wavelength bin.

Unique selling points:
- Includes a built-in simulation used to provide realistic uncertainties
- Easily applicable to different instruments and datasets
- Ideal for future large spectroscopic surveys such as the WEAVE-LOFAR survey ([Smith et al. 2016](https://ui.adsabs.harvard.edu/abs/2016sf2a.conf..271S/abstract))


If you make use of `SpecStacker` in your research, please include the following citation:

    @ARTICLE{2024MNRAS.528.4547A,
       author = {{Arnaudova}, M.~I. and {Smith}, D.~J.~B. and {Hardcastle}, M.~J. and {Das}, S. and {Drake}, A. and {Duncan}, K. and {G{\"u}rkan}, G. and {Magliocchetti}, M. and {Morabito}, L.~K. and {Petley}, J.~W. and {Shenoy}, S. and {Tasse}, C.},
        title = "{Exploring the radio loudness of SDSS quasars with spectral stacking}",
      journal = {\mnras},
     keywords = {techniques: spectroscopic, galaxies: active, quasars: general, radio continuum: galaxies, Astrophysics - Astrophysics of Galaxies},
         year = 2024,
        month = mar,
       volume = {528},
       number = {3},
        pages = {4547-4567},
          doi = {10.1093/mnras/stae233},
    archivePrefix = {arXiv},
       eprint = {2401.08774},
    primaryClass = {astro-ph.GA},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2024MNRAS.528.4547A},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }

