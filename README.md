# Whole-field-density-measurements-by-Digital-Image-Correlation

This repository contains the software used in "Whole-field density measurements by Digital Image Correlation", Alexander M. van Oers, Roeland de Kat, Leo R.M. Maas, 23 May 2023, PREPRINT (Version 1) available at Research Square [https://doi.org/10.21203/rs.3.rs-2947633/v1]. 

Abstract:
A new optical method is applied for the quantitative measurement of the density field of two-dimensional, stratified or homogeneous, transparent fluids in a laboratory set-up. A crucial new step is the angle at which we place a camera to view the experimental set-up. This step is motivated by a fallacy observed when applying our method in a classical configuration, in which the camera's optical axis is perpendicular to the flat surface of a fluid container. Application of this method is illustrated by the optical determination of static density fields of linearly and nonlinearly stratified fluids, as well as of multi-layered fluids. The method is validated by comparing with density profiles obtained from probe measurements of conductivity and temperature. Our method yields similar density and buoyancy profiles as the probe while also providing a whole field measurement, without disturbing the fluid and allowing the determination of dynamical density fields.

The software performs three main steps: 1) digital image correlation (dic), 2) calibration and (3) calculation of the index of refraction. The main file is DisplayImage.cpp. 

Expected usage of the software is "./DIC <Image_Path> SplineDegree SubsetLength GridLength ShapeFunction PropagationFunction OrderingImages xStart xEnd yStart yEnd NumberOfThreads MaxPixelYVertical Tolerance MinCorrCoeffIG BlurSize directionsToInclude nref DICNeeded CalibrationNeeded CalculateRefractionIndex", where <Image_Path> is the path to the folder where the images are located; SplineDegree, SubsetLength and GridLength are parameters of the dic step; ShapeFunction determines the order of the shape functions with 0 = rigid, 1 = affine, 2 = mixed terms, 3 = second order terms; PropagationFunction is a flag to use a propgationfunction to transfer the solution of the dic step of one grid point to another grid point as initial conditon. OrderingImages deterimes which image is the template and which image is the deformed image in the dic step. xStart, xEnd, yStart and yEnd indicate the size of the template image. NumberOfThreads is the number of 
	

As input two images are required, between which the dic will be performed. These images must be in the format *.tif or *.csv (in a folder called Averaged).


Tested on: Fedora Linux
Required Software: make, cmake, gcc, g++, openCV
