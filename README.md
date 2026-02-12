# WebSky Line Intensity Mapping

This repository contains the latest, stable version of the line intensity mapping (LIM) implementation of WebSky. These codes are based on the `limlam_mocker` package (originally written by George Stein, and later modified significantly by Dongwoo Chung, Patrick Breysse, Clara Chung and Patrick Horlaville). This implementation was forked from my branch of Patrick's CITA_LIM repo and streamlined to include just critical components. Currently, this code is being developed by Nate Carlson and Doga Tolgay to produce FIRE-informed CO mocks.

## Summary of code components

1. `lim.py` is the primary interface for interacting with the package. A lim model is loaded in python by importing `from lim import lim` and creating an object, e.g. `model = lim( model_params='COMAP_Fid', doSim=True )`
1. `mass_luminosity.py` contains the response functions to halo catalogues
1. ...
