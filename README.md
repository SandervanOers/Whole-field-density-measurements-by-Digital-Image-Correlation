# Whole-field-density-measurements-by-Digital-Image-Correlation

This repository contains the software used in "Whole-field density measurements by Digital Image Correlation", Alexander M. van Oers, Roeland de Kat, Leo R.M. Maas, 23 May 2023, PREPRINT (Version 1) available at Research Square [https://doi.org/10.21203/rs.3.rs-2947633/v1]. 

Abstract:
A new optical method is applied for the quantitative measurement of the density field of two-dimensional, stratified or homogeneous, transparent fluids in a laboratory set-up. A crucial new step is the angle at which we place a camera to view the experimental set-up. This step is motivated by a fallacy observed when applying our method in a classical configuration, in which the camera's optical axis is perpendicular to the flat surface of a fluid container. Application of this method is illustrated by the optical determination of static density fields of linearly and nonlinearly stratified fluids, as well as of multi-layered fluids. The method is validated by comparing with density profiles obtained from probe measurements of conductivity and temperature. Our method yields similar density and buoyancy profiles as the probe while also providing a whole field measurement, without disturbing the fluid and allowing the determination of dynamical density fields.

The software performs three main steps: 1) digital image correlation (dic), 2) calibration and (3) calculation of the index of refraction. The main file is DisplayImage.cpp. 

The expected usage of the software is explained in ExpectedUsage.md.

Tested on: Fedora Linux
Required Software: make, cmake, gcc, g++, openCV

Installation: Download repository to folder, "cd folder", "cmake .", "make ./DIC"
