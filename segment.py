# -*- coding: utf-8 -*-
"""
Creates class of functions required to perform manually seeded lesion segmentation

Created on Thu May 08 15:03:11 2014

@ author (C) Cristina Gallego, University of Toronto, 2014
"""

import os, os.path
import sys
import string
from sys import argv, stderr, exit
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import numpy as np

from display import *
import pylab

#!/usr/bin/env python
class Segment(object):
    """
    USAGE:
    =============
    loadSegment = Segment()
    """
    def __init__(self): 
        """ initialize visualization with standard vtk actors, renders, windowsn, interactors """           
        # use cell picker for interacting with the image orthogonal views.
        self.boxWidget = []
        self.boundsPlane_presel = []
        # Create only 1 display
        self.loadDisplay = Display()
        
    def __call__(self):       
        """ Turn Class into a callable object """
        Segment()
        
    def SelectPolygons(self, widget, event_string):
        # As can be seen the callback takes two arguments.  The first being the object that generates the event 
        # and the second argument the event name (which is a string).        
        self.planes = vtk.vtkPlanes()
        self.boxWidget.GetPlanes(self.planes)
        self.boundsPlane_presel = self.planes.GetPoints().GetBounds()
        
        return
    
    def ensureInSegment(self, image, lesionMesh, pathSegment, nameSegment, image_pos_pat, image_ori_pat):    
        
        # Proceed to build reference frame for display objects based on DICOM coords   
        [transformed_image, transform_cube] = self.loadDisplay.dicomTransform(image, image_pos_pat, image_ori_pat)
        
        dataToStencil = vtk.vtkPolyDataToImageStencil()
        dataToStencil.SetInput(lesionMesh)
        dataToStencil.SetOutputOrigin(transformed_image.GetOrigin())
        print transformed_image.GetOrigin()
        dataToStencil.SetOutputSpacing(transformed_image.GetSpacing())
        print transformed_image.GetSpacing()
        dataToStencil.SetOutputWholeExtent(transformed_image.GetExtent())
        dataToStencil.Update()
        
        stencil = vtk.vtkImageStencil()
        stencil.SetInput(transformed_image)
        stencil.SetStencil(dataToStencil.GetOutput())
        stencil.ReverseStencilOff()
        stencil.SetBackgroundValue(0.0)
        stencil.Update()
        
        newSegment = vtk.vtkMetaImageWriter()
        newSegment.SetFileName(pathSegment+'/'+nameSegment+'.mhd')
        newSegment.SetInput(stencil.GetOutput())
        newSegment.Write()

        thresh = vtk.vtkImageThreshold()
        thresh.SetInput(stencil.GetOutput())
        thresh.ThresholdByUpper(1)
        thresh.SetInValue(255)
        thresh.SetOutValue(0)
        thresh.Update()
                  
        contouriso = vtk.vtkMarchingCubes()
        contouriso.SetInput(thresh.GetOutput())
        contouriso.SetValue(0,125)
        contouriso.ComputeScalarsOn()
        contouriso.Update()
        
        # Recalculate num_voxels and vol_lesion on VOI
        nvoxels = contouriso.GetOutput().GetNumberOfCells()
        npoints = contouriso.GetOutput().GetNumberOfPoints()
        print "Number of points: %d" % npoints 
        print "Number of cells: %d" % nvoxels 
        
        return contouriso.GetOutput()
        
        
    def saveSegmentation(self, lesionID_path, lesionMesh, lesionfilename):
        """
        Saves new segmentation to file
        ARGUMENTS:
        =============
        lesionID_path: (str)        path to lesion segmentation
        lesionMesh (vtkPolyData)      3D lesion segmentation as a vtkPolyData object       
        """ 
        # need to locate and read Lesion Seg
        if lesionfilename:
            VOIlesion = lesionfilename    
        else:
            VOIlesion = 'VOIlesion_selected.vtk'
            
        lesion3D_writer = vtk.vtkPolyDataWriter()
        lesion3D_writer.SetInput(lesionMesh)
        lesion3D_writer.SetFileName( lesionID_path+os.sep+VOIlesion )
        lesion3D_writer.Update()
              
        return 
        
        
    def segmentFromSeeds(self, images,  image_pos_pat, image_ori_pat, seeds, iren, xplane, yplane, zplane):
        """ segmentFromSeeds: Extracts VOI from seeds
        
        INPUTS:
        =======        
        images: (vtkImageData)   list of Input image to Transform
        seeds: vtkPoints    list of seeded coordinates from picker
        OUTPUTS:
        =======
        transformed_image (vtkImageData)    Transformed imaged mapped to dicom coords frame
        transform (vtkTransform)            Transform used
        
        """
        for postS in range(1,len(images)):
            print "\n Segmenting image post no : %s " % str(postS)
            
            subimage = self.loadDisplay.subImage(images, postS)            
            
            # Proceed to build reference frame for display objects based on DICOM coords   
            [transformed_image, transform_cube] = self.loadDisplay.dicomTransform(subimage, image_pos_pat, image_ori_pat)
            
            # Calculate the center of the volume
            transformed_image.UpdateInformation() 
                    
            print "\nBoxwidget placed..."
            #################################################################
            # The box widget observes the events invoked by the render window
            # interactor.  These events come from user interaction in the render
            # window.
            # Place the interactor initially. The output of the reader is used to
            # place the box widget.
            self.boxWidget = vtk.vtkBoxWidget()
            self.boxWidget.SetInteractor(iren)
            self.boxWidget.SetPlaceFactor(1)
            self.boxWidget.SetInput(transformed_image)
            if( self.boundsPlane_presel != []):
                self.boxWidget.PlaceWidget( self.boundsPlane_presel )
            
            # Construct a bounding box around the seeds  
            init_seedsBounds = [0,0,0,0,0,0]
            seeds.GetBounds( init_seedsBounds )
            
            if postS == 1:        
                # polygonal data --> image stencil:
                # create a simple box VOI mask shape using previously found boundsPlane_preselected
                VOIStencil = vtk.vtkROIStencilSource()
                VOIStencil.SetShapeToBox()
                VOIStencil.SetBounds( init_seedsBounds )    
                VOIStencil.SetInformationInput(transformed_image)
                VOIStencil.Update()
            
                # cut the corresponding VOI region and set the background:
                extractVOI_imgstenc = vtk.vtkImageStencil()
                extractVOI_imgstenc.SetInput(transformed_image)
                extractVOI_imgstenc.SetStencil(VOIStencil.GetOutput())
                extractVOI_imgstenc.ReverseStencilOff()
                extractVOI_imgstenc.SetBackgroundValue(0.0)
                extractVOI_imgstenc.Update()
                
                allSeededIm = vtk.vtkImageData()
                allSeededIm.DeepCopy( extractVOI_imgstenc.GetOutput() )
        
                # Add some bounding box radius
                self.boxWidget.PlaceWidget( init_seedsBounds )
                self.boundsPlane_presel = init_seedsBounds
                print "seeds.GetBounds"
                print init_seedsBounds
                        
                self.boxWidget.AddObserver("InteractionEvent", self.SelectPolygons)
                self.boxWidget.On()
            
                # turn off planes
                xplane.Off()
                yplane.Off()
                iren.Start()
                self.boxWidget.Off()
            
            # polygonal data --> image stencil:
            print "\n Create vtkPolyDataToImageStencil with bounds:"
            print self.boundsPlane_presel
            
            # create a simple box VOI mask shape using previously found boundsPlane_preselected
            VOIStencil = vtk.vtkROIStencilSource()
            VOIStencil.SetShapeToBox()
            VOIStencil.SetBounds( self.boundsPlane_presel )    
            VOIStencil.SetInformationInput(transformed_image)
            VOIStencil.Update()
                                    
            # cut the corresponding VOI region and set the background:
            extractVOI_imgstenc = vtk.vtkImageStencil()
            extractVOI_imgstenc.SetInput(transformed_image)
            extractVOI_imgstenc.SetStencil(VOIStencil.GetOutput())
            extractVOI_imgstenc.ReverseStencilOff()
            extractVOI_imgstenc.SetBackgroundValue(0.0)
            extractVOI_imgstenc.Update()
                
            # add subsecuent VOI stencils
            addsegROI = vtk.vtkImageMathematics()
            addsegROI.SetInput(0, allSeededIm)
            addsegROI.SetInput(1, extractVOI_imgstenc.GetOutput())
            addsegROI.SetOperationToAdd()
            addsegROI.Update()
      
            # turn off the box
            allSeededIm = addsegROI.GetOutput()

        avegSeededIm = vtk.vtkImageMathematics()
        avegSeededIm.SetInput(0, allSeededIm)
        avegSeededIm.SetOperationToMultiplyByK()
        avegSeededIm.SetConstantK( 0.2 )        
        avegSeededIm.Update()
        
        # take out average image
        finalSeedIm = avegSeededIm.GetOutput()
        xplane.SetInput( finalSeedIm )
        yplane.SetInput( finalSeedIm )
        zplane.SetInput( finalSeedIm )            
        
        image_scalar_range = finalSeedIm.GetScalarRange() 
        lThre = image_scalar_range[0]
        uThre = image_scalar_range[1]
        print "\n Image Scalar Range:"
        print image_scalar_range[0], image_scalar_range[1]
        print "Uper thresholding by"
        print uThre*0.25
        
        ## Display histogram 
        scalars = finalSeedIm.GetPointData().GetScalars()
        np_scalars = vtk_to_numpy(scalars)        
        pylab.figure()
        pylab.hist(np_scalars.flatten(), histtype='bar')
        pylab.show()
                
        # Now performed threshold on initialized VOI        
        # vtkImageThresholdConnectivity will perform a flood fill on an image, given upper and lower pixel intensity 
        # thresholds. It works similarly to vtkImageThreshold, but also allows the user to set seed points to limit
        # the threshold operation to contiguous regions of the image. The filled region, or the "inside", will be passed 
        # through to the output by default, while the "outside" will be replaced with zeros. The scalar type of the output is the same as the input.
        thresh_sub = vtk.vtkImageThresholdConnectivity() 
        thresh_sub.SetSeedPoints(seeds)
        thresh_sub.SetNeighborhoodRadius(3, 3, 2) #The radius of the neighborhood that must be within the threshold values in order for the voxel to be included in the mask. The default radius is zero (one single voxel). The radius is measured in voxels
        thresh_sub.SetNeighborhoodFraction(0.10) #The fraction of the neighborhood that must be within the thresholds. The default value is 0.5.
        thresh_sub.ThresholdBetween(0.25*uThre, uThre); 
        thresh_sub.SetInput( finalSeedIm )
        thresh_sub.Update()
        
        xplane.SetInput( thresh_sub.GetOutput() )
        yplane.SetInput( thresh_sub.GetOutput()  )
        zplane.SetInput( thresh_sub.GetOutput()  )            
        iren.Start()
        
        # Convert VOIlesion into polygonal struct
        VOIlesion_poly = vtk.vtkMarchingCubes() 
        VOIlesion_poly.SetValue(0,255)
        VOIlesion_poly.SetInput(thresh_sub.GetOutput())
        VOIlesion_poly.ComputeNormalsOff()
        VOIlesion_poly.Update()
        
        # Recalculate num_voxels and vol_lesion on VOI
        nvoxels = VOIlesion_poly.GetOutput().GetNumberOfCells()
        print "\n Number of Voxels: %d" % nvoxels # After the filter has executed, use GetNumberOfVoxels() to find out how many voxels were filled.

        return VOIlesion_poly.GetOutput()
