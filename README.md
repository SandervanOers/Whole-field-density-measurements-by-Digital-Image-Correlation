# Whole-field-density-measurements-by-Digital-Image-Correlation

This repository contains the software used in "Whole-field density measurements by Digital Image Correlation", Alexander M. van Oers, Roeland de Kat, Leo R.M. Maas, 23 May 2023, PREPRINT (Version 1) available at Research Square [https://doi.org/10.21203/rs.3.rs-2947633/v1]. 

Abstract:
A novel application of Synthetic Schlieren in a laboratory set-up yields a quantitative measurement of the density field of two-dimensional, stratified or homogeneous, transparent fluids in a laboratory set-up using a single camera. This application obtains local values of the density without the need for tomographic reconstruction algorithms that require images taken from different directions through the fluid nor does the application require regularization. This is achieved by placing the camera at a large oblique angle with respect to the experimental set-up. This step is motivated by a fallacy observed when applying ray tracing in a classical configuration, in which the camera's optical axis is perpendicular to the flat surface of a fluid container. The application is illustrated by the optical determination of static density fields of linearly and nonlinearly stratified fluids, as well as of multi-layered fluids. The application is validated by comparing with density profiles obtained from probe measurements of conductivity and temperature. Our application yields similar density and density gradient profiles as the probe while also providing a whole field measurement without disturbing the fluid,	 and allowing the determination of dynamical density fields.

The software performs three main steps: 1) digital image correlation (dic), 2) calibration and (3) calculation of the index of refraction. The main file is DisplayImage.cpp. 

The expected usage of the software is explained in ExpectedUsage.md.

Tested on: Fedora Linux
Required Software: make, cmake, gcc, g++, openCV

Installation: Download repository to folder, "cd folder", "cmake .", "make ./DIC"
