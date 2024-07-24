# SpecStacker: The Spectral Stacking Code
This is a new rest-frame spectral stacking code.

The code takes an array of observed-frame wavelength, flux density, flux uncertainties and the corresponding redshifts as input parameters, shifts all spectra to the rest-frame, and resamples them onto a common wavelegth grid based on the redshift distribution and an input

and outputs the best-fit parameters and quality-checking plots to the paths specified by the user.

If you make use of `SpecStacker` in your research please include the following citation:

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

