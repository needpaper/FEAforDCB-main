# FEAforDCB.py
# ---------------------------------------
# This program is a graphical user interface (GUI) for Finite Element Analysis (FEA) tasks.
# Users can configure model parameters, select analysis types, and generate corresponding models,
# or perform static analysis, nonlinear buckling analysis, and frequency analysis.
#
# Contents:
# 1. GUI Initialization and Layout: Defines the window, buttons, input fields, and other interface elements.
# 2. Functional Module Integration: Includes model creation, material setting, mesh generation,
#    and the invocation of various analysis modules.
# 3. Event Handling: Binds GUI events (e.g., button clicks) to corresponding functional modules.
#
#
# Writen by Wang Haobo .  HIT .Harbin .China
# Date: October 2024
#
#
from abaqus import *
from abaqusConstants import *
from part import *
from material import *
import section
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
from visualization import *
from connectorBehavior import *
from abaqus import session
from odbAccess import *
from odbSection import *
import numpy as np
import math
import os
import csv


def runFEAforDCB(model_lib, Design_space, length, radius, angle ,materialdensity,materialtable,
              materialsection ,meshsize,anslysislib,Staticnlgeom,Bendingdirection,StaticCPUnum,Staticsaveresult,
              Bucklingdirection,Bucklingmode,defectcoefficient,saveresult,NonLinearCPUnum,Eigenreq,RefBC,Inertiamass,freqCPUnum):

   

    model = mdb.models['Model-1']

    def Modeling(model_lib, Design_space, length, radius, angle,materialdensity, materialtable, materialsection,meshsize):
        if model_lib == 'CTLT':
            partName = 'CTLT'
            l1 = length
            r1 = radius
            a1 = angle
            r11 = r1 * cos(angle * pi / 180.0) / (1 - cos(angle * pi / 180.0))
            w1 = (Design_space - angle * r1 * pi / 90 - angle * r11 * pi / 90) / 2
            ##Part##
            model.ConstrainedSketch(name='profile', sheetSize=400.0)
            model.sketches['profile'].ArcByCenterEnds(
                center=(0.0, 0.0),
                direction=CLOCKWISE,
                point1=(-r1 * sin(a1 * pi / 180.0), r1 * cos(a1 * pi / 180.0)),
                point2=(r1 * sin(a1 * pi / 180.0), r1 * cos(a1 * pi / 180.0)))
            model.sketches['profile'].ArcByCenterEnds(
                center=(r1 * sin(a1 * pi / 180.0) + r11 * sin(a1 * pi / 180.0), r11),
                direction=COUNTERCLOCKWISE,
                point1=(r1 * sin(a1 * pi / 180.0), r1 * cos(a1 * pi / 180.0)),
                point2=(r1 * sin(a1 * pi / 180.0) + r11 * sin(a1 * pi / 180.0), 0.0))
            model.sketches['profile'].ArcByCenterEnds(
                center=(-r1 * sin(a1 * pi / 180.0) - r11 * sin(a1 * pi / 180.0), r11),
                direction=CLOCKWISE,
                point1=(-r1 * sin(a1 * pi / 180.0), r1 * cos(a1 * pi / 180.0)),
                point2=(-r1 * sin(a1 * pi / 180.0) - r11 * sin(a1 * pi / 180.0), 0.0))
            model.sketches['profile'].Line(
                point1=(r1 * sin(a1 * pi / 180.0) + r11 * sin(a1 * pi / 180.0), 0.0),
                point2=(r1 * sin(a1 * pi / 180.0) + r11 * sin(a1 * pi / 180.0) + w1, 0.0))
            model.sketches['profile'].Line(
                point1=(-r1 * sin(a1 * pi / 180.0) - r11 * sin(a1 * pi / 180.0), 0.0),
                point2=(-r1 * sin(a1 * pi / 180.0) - r11 * sin(a1 * pi / 180.0) - w1, 0.0))
            model.sketches['profile'].ArcByCenterEnds(
                center=(0.0, 0.0),
                direction=CLOCKWISE,
                point1=(r1 * sin(a1 * pi / 180.0), -r1 * cos(a1 * pi / 180.0)),
                point2=(-r1 * sin(a1 * pi / 180.0), -r1 * cos(a1 * pi / 180.0)))
            model.sketches['profile'].ArcByCenterEnds(
                center=(r1 * sin(a1 * pi / 180.0) + r11 * sin(a1 * pi / 180.0), -r11),
                direction=CLOCKWISE,
                point1=(r1 * sin(a1 * pi / 180.0), -r1 * cos(a1 * pi / 180.0)),
                point2=(r1 * sin(a1 * pi / 180.0) + r11 * sin(a1 * pi / 180.0), 0.0))
            model.sketches['profile'].ArcByCenterEnds(
                center=(-r1 * sin(a1 * pi / 180.0) - r11 * sin(a1 * pi / 180.0), -r11),
                direction=COUNTERCLOCKWISE,
                point1=(-r1 * sin(a1 * pi / 180.0), -r1 * cos(a1 * pi / 180.0)),
                point2=(-r1 * sin(a1 * pi / 180.0) - r11 * sin(a1 * pi / 180.0), 0.0))
            model.Part(dimensionality=THREE_D, name=partName, type=
                DEFORMABLE_BODY)
            model.parts[partName].BaseShellExtrude(depth=l1, sketch=
                model.sketches['profile'])
            del model.sketches['profile']

            print('A CTLT part is created')

            print('Setting material for CTLT part')
            ##Material##
            ##set ctlt_single and ctlt_double
            model.parts[partName].Set(faces=
            model.parts[partName].faces.findAt((
                (r1 * sin(a1 * pi / 180.0) + r11 * sin(a1 * pi / 180.0) - r11 * sin(a1 * pi / 360.0),
                 r11 - r11 * cos(a1 * pi / 360.0),
                 l1 / 2),), ((0.0, r1, l1 / 2),),
                ((-(r1 * sin(a1 * pi / 180.0) + r11 * sin(a1 * pi / 180.0) - r11 * sin(a1 * pi / 360.0)),
                  r11 - r11 * cos(a1 * pi / 360.0), l1 / 2),),
                ((r1 * sin(a1 * pi / 180.0) + r11 * sin(a1 * pi / 180.0) - r11 * sin(a1 * pi / 360.0),
                  -(r11 - r11 * cos(a1 * pi / 360.0)), l1 / 2),), ((0.0, -r1, l1 / 2),), (
                    (-(r1 * sin(a1 * pi / 180.0) + r11 * sin(a1 * pi / 180.0) - r11 * sin(a1 * pi / 360.0)),
                     -(r11 - r11 * cos(a1 * pi / 360.0)),
                     l1 / 2),), ), name='ctlt_single')

            model.parts[partName].Set(faces=
                                      model.parts[
                                          partName].faces.findAt(((r1 * sin(
                                          a1 * pi / 180.0) + r11 * sin(a1 * pi / 180.0) + w1,
                                                                   0.0, 0.0),), ((r1 * sin(
                                          a1 * pi / 180.0) + r11 * sin(a1 * pi / 180.0) + w1,
                                                                                  0.0, l1),), ),
                                      name='ctlt_doub1')
            model.parts[partName].Set(faces=
                                      model.parts[
                                          partName].faces.findAt(((-r1 * sin(
                                          a1 * pi / 180.0) - r11 * sin(a1 * pi / 180.0) - w1,
                                                                   0.0, 0.0),), ((-r1 * sin(
                                          a1 * pi / 180.0) - r11 * sin(a1 * pi / 180.0) - w1,
                                                                                  0.0, l1),), ),
                                      name='ctlt_doub2')

            ##material
            model.Material(name='Composite')
            model.materials['Composite'].Elastic(table=materialtable,
                                                 type=LAMINA)
            model.materials['Composite'].Density(table=((materialdensity,),))

            ##SectionAssignment

            layup_tuple = tuple(
                section.SectionLayer(
                    material='Composite',
                    thickness=float(layer[0]),
                    orientAngle=float(layer[1])
                ) for layer in materialsection
            )
            model.CompositeShellSection(name='Section-single',
                                        idealization=NO_IDEALIZATION,
                                        integrationRule=SIMPSON,
                                        layup=layup_tuple,
                                        poissonDefinition=DEFAULT, preIntegrate=OFF,
                                        symmetric=False, temperature=GRADIENT,
                                        thicknessModulus=None, thicknessType=UNIFORM,
                                        useDensity=OFF)

            model.parts[partName].SectionAssignment(
                offset=0.0,
                offsetField='',
                offsetType=MIDDLE_SURFACE,
                region=model.parts[partName].sets['ctlt_single'],
                sectionName='Section-single',
                thicknessAssignment=FROM_SECTION)

            symmetric_layup = layup_tuple + tuple(
                section.SectionLayer(
                    material='Composite',
                    thickness=layer.thickness,
                    orientAngle=layer.orientAngle
                ) for layer in reversed(layup_tuple)
            )
            model.CompositeShellSection(
                name='Section-double',
                idealization=NO_IDEALIZATION,
                integrationRule=SIMPSON,
                layup=symmetric_layup,
                poissonDefinition=DEFAULT, preIntegrate=OFF,
                symmetric=True, temperature=GRADIENT,
                thicknessModulus=None, thicknessType=UNIFORM,
                useDensity=OFF)

            model.parts[partName].SectionAssignment(
                offset=0.0,
                offsetField='',
                offsetType=MIDDLE_SURFACE,
                region=model.parts[partName].sets['ctlt_doub1'],
                sectionName='Section-double',
                thicknessAssignment=FROM_SECTION)
            model.parts[partName].SectionAssignment(
                offset=0.0,
                offsetField='',
                offsetType=MIDDLE_SURFACE,
                region=model.parts[partName].sets['ctlt_doub2'],
                sectionName='Section-double',
                thicknessAssignment=FROM_SECTION)

            model.parts[partName].MaterialOrientation(additionalRotationType=
                                                      ROTATION_NONE, axis=AXIS_2, fieldName='',
                                                      localCsys=None, orientationType=
                                                      GLOBAL, region=
                                                      model.parts[
                                                          partName].sets['ctlt_single'])
            model.parts[partName].MaterialOrientation(additionalRotationType=
                                                      ROTATION_NONE, axis=AXIS_2, fieldName='',
                                                      localCsys=None, orientationType=
                                                      GLOBAL, region=
                                                      model.parts[
                                                          partName].sets['ctlt_doub1'])
            model.parts[partName].MaterialOrientation(additionalRotationType=
                                                      ROTATION_NONE, axis=AXIS_2, fieldName='',
                                                      localCsys=None, orientationType=
                                                      GLOBAL, region=
                                                      model.parts[
                                                          partName].sets['ctlt_doub2'])

            print('Assembling and Meshing for CTLT part')
            ##Assmebly##
            model.rootAssembly.DatumCsysByDefault(CARTESIAN)
            model.rootAssembly.Instance(dependent=ON, name='ctlt-1', part=
                 model.parts[partName])

            ##creat referencePoints
            model.rootAssembly.ReferencePoint(point=(0.0, 0.0, 0.0))
            model.rootAssembly.ReferencePoint(point=(0.0, 0.0, l1))

            ##creat section
            model.rootAssembly.Set(edges=
                model.rootAssembly.instances['ctlt-1'].edges.findAt(((
                                                                     -(r1 * sin(
                                                                         a1 * pi / 180.0) + r11 * sin(
                                                                         a1 * pi / 180.0) + w1 / 2),
                                                                     0.0, 0.0),),
                ((-(r1 * sin(a1 * pi / 180.0) + r11 * sin(a1 * pi / 180.0) - r11 * sin(a1 * pi / 360.0)),
                  r11 - r11 * cos(a1 * pi / 360.0), 0.0),),
                ((0.0, r1, 0.0),),
                ((r1 * sin(a1 * pi / 180.0) + r11 * sin(a1 * pi / 180.0) - r11 * sin(a1 * pi / 360.0),
                  r11 - r11 * cos(a1 * pi / 360.0), 0.0),),
                ((r1 * sin(a1 * pi / 180.0) + r11 * sin(a1 * pi / 180.0) + w1 / 2, 0.0, 0.0),),
                ((-(r1 * sin(a1 * pi / 180.0) + r11 * sin(a1 * pi / 180.0) - r11 * sin(a1 * pi / 360.0)),
                  -(r11 - r11 * cos(a1 * pi / 360.0)), 0.0),),
                ((0.0, -r1, 0.0),),
                ((r1 * sin(a1 * pi / 180.0) + r11 * sin(a1 * pi / 180.0) - r11 * sin(a1 * pi / 360.0),
                  -(r11 - r11 * cos(a1 * pi / 360.0)), 0.0),)), name='ctlt_section1')
            model.rootAssembly.Set(edges=
               model.rootAssembly.instances['ctlt-1'].edges.findAt(((
                                                                     -(r1 * sin(
                                                                         a1 * pi / 180.0) + r11 * sin(
                                                                         a1 * pi / 180.0) + w1 / 2),
                                                                     0.0, l1),),
                ((-(r1 * sin(a1 * pi / 180.0) + r11 * sin(a1 * pi / 180.0) - r11 * sin(a1 * pi / 360.0)),
                  r11 - r11 * cos(a1 * pi / 360.0), l1),),
                ((0.0, r1, l1),),
                ((r1 * sin(a1 * pi / 180.0) + r11 * sin(a1 * pi / 180.0) - r11 * sin(a1 * pi / 360.0),
                  r11 - r11 * cos(a1 * pi / 360.0), l1),),
                ((r1 * sin(a1 * pi / 180.0) + r11 * sin(a1 * pi / 180.0) + w1 / 2, 0.0, l1),),
                ((-(r1 * sin(a1 * pi / 180.0) + r11 * sin(a1 * pi / 180.0) - r11 * sin(a1 * pi / 360.0)),
                  -(r11 - r11 * cos(a1 * pi / 360.0)), l1),),
                ((0.0, -r1, l1),),
                ((r1 * sin(a1 * pi / 180.0) + r11 * sin(a1 * pi / 180.0) - r11 * sin(a1 * pi / 360.0),
                  -(r11 - r11 * cos(a1 * pi / 360.0)), l1),)), name='ctlt_section2')

            ##Interaction##
            ##create coupling RP_with_ctltsection
            model.Coupling(name='RP-1_with_ctltsection1', controlPoint=
                Region(referencePoints=(model.rootAssembly.referencePoints[4],)),
                           surface=model.rootAssembly.sets[
                               'ctlt_section1'],
                           influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC,
                           localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)
            model.Coupling(name='RP-2_with_ctltsection2', controlPoint=
                Region(referencePoints=(model.rootAssembly.referencePoints[5],)),
                           surface=model.rootAssembly.sets[
                               'ctlt_section2'],
                           influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC,
                           localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)
            ##mesh##
            model.parts[partName].seedPart(deviationFactor=0.1, minSizeFactor=
            0.1, size=meshsize)
            model.parts[partName].generateMesh()
            model.parts[partName].setElementType(elemTypes=(
                ElemType(elemCode=S4R, elemLibrary=EXPLICIT, secondOrderAccuracy=OFF, hourglassControl=STIFFNESS),
                ElemType(elemCode=S3R, elemLibrary=EXPLICIT)), regions=(model.parts[partName].faces,))
            ##set node
            model.rootAssembly.Set(referencePoints=(model.rootAssembly.referencePoints[4],), name='Ref1')
            model.rootAssembly.Set(referencePoints=(model.rootAssembly.referencePoints[5],), name='Ref2')


        if model_lib == 'TRAC':
            partName = 'TRAC'
            l2 = length
            w2 = Design_space - angle * radius * pi / 180
            r2 = radius
            a2 = angle
            ##Part##
            model.ConstrainedSketch(name='profile', sheetSize=400.0)
            model.sketches['profile'].Line(
                point1=(0.0, 0.0),
                point2=(w2, 0.0))
            model.sketches['profile'].ArcByCenterEnds(
                center=(0.0, r2),
                direction=CLOCKWISE,
                point1=(0.0, 0.0),
                point2=(-r2 * sin(a2 * pi / 180.0), r2 - r2 * cos(a2 * pi / 180.0)))
            model.sketches['profile'].ArcByCenterEnds(
                center=(0.0, -r2),
                direction=CLOCKWISE,
                point1=(-r2 * sin(a2 * pi / 180.0), r2 * cos(a2 * pi / 180.0) - r2),
                point2=(0.0, 0.0))

            model.Part(dimensionality=THREE_D, name=partName, type=
                DEFORMABLE_BODY)
            model.parts[partName].BaseShellExtrude(depth=l2, sketch=
                model.sketches['profile'])
            del model.sketches['profile']

            print('A TRAC part is created')

            print('Setting material for TRAC part')

            ##Material##
            ##set TRAC_single and TRAC_double
            model.parts[partName].Set(faces=
            model.parts[partName].faces.findAt(
                ((0.5 * w2, 0.0, l2 / 2),), ), name='TRAC_double')
            model.parts[partName].Set(faces=
                                      model.parts[partName].faces.findAt(((-r2 * sin(
                                          a2 * pi / 180.0), r2 - r2 * cos(a2 * pi / 180.0),
                                                                           l2 / 2),), ),
                                      name='TRAC_single1')
            model.parts[partName].Set(faces=
                                      model.parts[partName].faces.findAt(((-r2 * sin(
                                          a2 * pi / 180.0), r2 * cos(a2 * pi / 180.0) - r2,
                                                                           l2 / 2),), ),
                                      name='TRAC_single2')

            ##material
            model.Material(name='Composite')
            model.materials['Composite'].Elastic(table=materialtable,
                                                 type=LAMINA)
            model.materials['Composite'].Density(table=((materialdensity,),))

            ##SectionAssignment
            layup_tuple = tuple(
                section.SectionLayer(
                    material='Composite',
                    thickness=float(layer[0]),
                    orientAngle=float(layer[1])
                ) for layer in materialsection
            )
            model.CompositeShellSection(name='Section-single',
                                        idealization=NO_IDEALIZATION,
                                        integrationRule=SIMPSON,
                                        layup=layup_tuple,
                                        poissonDefinition=DEFAULT, preIntegrate=OFF,
                                        symmetric=False, temperature=GRADIENT,
                                        thicknessModulus=None, thicknessType=UNIFORM,
                                        useDensity=OFF)

            model.parts[partName].SectionAssignment(
                offset=0.0,
                offsetField='',
                offsetType=MIDDLE_SURFACE,
                region=model.parts[partName].sets['TRAC_single1'],
                sectionName='Section-single',
                thicknessAssignment=FROM_SECTION)
            model.parts[partName].SectionAssignment(
                offset=0.0,
                offsetField='',
                offsetType=MIDDLE_SURFACE,
                region=model.parts[partName].sets['TRAC_single2'],
                sectionName='Section-single',
                thicknessAssignment=FROM_SECTION)
            symmetric_layup = layup_tuple + tuple(
                section.SectionLayer(
                    material='Composite',
                    thickness=layer.thickness,
                    orientAngle=layer.orientAngle
                ) for layer in reversed(layup_tuple)
            )
            model.CompositeShellSection(
                name='Section-double',
                idealization=NO_IDEALIZATION,
                integrationRule=SIMPSON,
                layup=symmetric_layup,
                poissonDefinition=DEFAULT, preIntegrate=OFF,
                symmetric=True, temperature=GRADIENT,
                thicknessModulus=None, thicknessType=UNIFORM,
                useDensity=OFF)

            model.parts[partName].SectionAssignment(
                offset=0.0,
                offsetField='',
                offsetType=MIDDLE_SURFACE,
                region=model.parts[partName].sets['TRAC_double'],
                sectionName='Section-double',
                thicknessAssignment=FROM_SECTION)

            model.parts[partName].MaterialOrientation(additionalRotationType=
                                                      ROTATION_NONE, axis=AXIS_2, fieldName='',
                                                      localCsys=None, orientationType=
                                                      GLOBAL,
                                                      region=
                                                      model.parts[partName].sets[
                                                          'TRAC_single1'])
            model.parts[partName].MaterialOrientation(additionalRotationType=
                                                      ROTATION_NONE, axis=AXIS_2, fieldName='',
                                                      localCsys=None, orientationType=
                                                      GLOBAL,
                                                      region=
                                                      model.parts[partName].sets[
                                                          'TRAC_single2'])
            model.parts[partName].MaterialOrientation(additionalRotationType=
                                                      ROTATION_NONE, axis=AXIS_2, fieldName='',
                                                      localCsys=None, orientationType=
                                                      GLOBAL,
                                                      region=
                                                      model.parts[partName].sets[
                                                          'TRAC_double'])

            print('Assembling and Meshing for TRAC part')

            ##Assmebly##
            model.rootAssembly.DatumCsysByDefault(CARTESIAN)
            model.rootAssembly.Instance(dependent=ON, name='TRAC-1', part=
            model.parts[partName])
            ##creat referencePoints
            model.rootAssembly.ReferencePoint(point=(0.0, 0.0, 0.0))
            model.rootAssembly.ReferencePoint(point=(0.0, 0.0, l2))

            ##creat section
            model.rootAssembly.Set(edges=
            model.rootAssembly.instances['TRAC-1'].edges.findAt(
                ((0.0, 0.0, 0.0),),
                ((-r2 * sin(a2 * pi / 360.0), r2 - r2 * cos(a2 * pi / 360.0),
                  0.0),),
                ((-r2 * sin(a2 * pi / 360.0), r2 * cos(a2 * pi / 360.0) - r2,
                  0.0),),
                ((w2 / 2, 0.0, 0.0),), ), name='TRAC_section1')
            model.rootAssembly.Set(edges=
            model.rootAssembly.instances['TRAC-1'].edges.findAt(
                ((0.0, 0.0, l2),),
                ((-r2 * sin(a2 * pi / 360.0), r2 - r2 * cos(a2 * pi / 360.0),
                  l2),),
                ((-r2 * sin(a2 * pi / 360.0), r2 * cos(a2 * pi / 360.0) - r2,
                  l2),),
                ((w2 / 2, 0.0, l2),), ), name='TRAC_section2')

            ##Interaction##
            ##create coupling RP_with_TRACsection
            model.Coupling(name='RP-1_with_TRACsection1', controlPoint=
                Region(referencePoints=(model.rootAssembly.referencePoints[4],)),
                           surface=model.rootAssembly.sets['TRAC_section1'],
                           influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC,
                           localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)
            model.Coupling(name='RP-2_with_TRACsection2', controlPoint=
                Region(referencePoints=(model.rootAssembly.referencePoints[5],)),
                           surface=model.rootAssembly.sets['TRAC_section2'],
                           influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC,
                           localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)

            ##mesh##
            model.parts[partName].seedPart(deviationFactor=0.1, minSizeFactor=
            0.1, size=meshsize)
            model.parts[partName].generateMesh()
            model.parts[partName].setElementType(elemTypes=(
                ElemType(elemCode=S4R, elemLibrary=EXPLICIT, secondOrderAccuracy=OFF, hourglassControl=STIFFNESS),
                ElemType(elemCode=S3R, elemLibrary=EXPLICIT)), regions=(model.parts[partName].faces,))
            ##set node
            model.rootAssembly.Set(referencePoints=(model.rootAssembly.referencePoints[4],), name='Ref1')
            model.rootAssembly.Set(referencePoints=(model.rootAssembly.referencePoints[5],), name='Ref2')


        if model_lib == 'CBOOM':
            partName = 'CBOOM'
            l3 = length
            a3 = angle
            r3 = Design_space / ((360 - a3) * pi / 180)

            ##Part##
            model.ConstrainedSketch(name='profile', sheetSize=400.0)
            model.sketches['profile'].ArcByCenterEnds(
                center=(0.0, 0.0),
                direction=CLOCKWISE,
                point1=(-r3 * sin(a3 * pi / 360.0), -r3 * cos(a3 * pi / 360.0)),
                point2=(r3 * sin(a3 * pi / 360.0), -r3 * cos(a3 * pi / 360.0)))

            model.Part(dimensionality=THREE_D, name=partName, type=
                DEFORMABLE_BODY)
            model.parts[partName].BaseShellExtrude(depth=l3, sketch=
                model.sketches['profile'])
            del model.sketches['profile']

            print('A CBOOM part is created')

            print('Setting material for CBOOM part')
            ##Material##
            model.parts[partName].Set(faces=
            model.parts[partName].faces.findAt(
                ((0.0, r3, 0.0),), ((0.0, r3, l3),), ), name='cboom_section')

            ##create material
            model.Material(name='Composite')
            model.materials['Composite'].Elastic(table=materialtable,
                                                 type=LAMINA)
            model.materials['Composite'].Density(table=((materialdensity,),))

            ##SectionAssignment
            layup_tuple = tuple(
                section.SectionLayer(
                    material='Composite',
                    thickness=float(layer[0]),
                    orientAngle=float(layer[1])
                ) for layer in materialsection
            )
            model.CompositeShellSection(name='Section-single',
                                        idealization=NO_IDEALIZATION,
                                        integrationRule=SIMPSON,
                                        layup=layup_tuple,
                                        poissonDefinition=DEFAULT, preIntegrate=OFF,
                                        symmetric=False, temperature=GRADIENT,
                                        thicknessModulus=None, thicknessType=UNIFORM,
                                        useDensity=OFF)

            model.parts[partName].SectionAssignment(
                offset=0.0,
                offsetField='',
                offsetType=MIDDLE_SURFACE,
                region=model.parts[partName].sets['cboom_section'],
                sectionName='Section-single',
                thicknessAssignment=FROM_SECTION)

            model.parts[partName].MaterialOrientation(additionalRotationType=
                                                      ROTATION_NONE, axis=AXIS_2, fieldName='',
                                                      localCsys=None, orientationType=
                                                      GLOBAL, region=
                                                      model.parts[partName].sets[
                                                          'cboom_section'])

            print('Assembling and Meshing for CBOOM part')

            ##Assmebly##
            model.rootAssembly.DatumCsysByDefault(CARTESIAN)
            model.rootAssembly.Instance(dependent=ON, name='cboom-1', part=
            model.parts[partName])

            ##creat referencePoints
            model.rootAssembly.ReferencePoint(point=(0.0, 0.0, 0.0))
            model.rootAssembly.ReferencePoint(point=(0.0, 0.0, l3))

            ##creat section
            model.rootAssembly.Set(edges=
            model.rootAssembly.instances['cboom-1'].edges.findAt(
                ((0.0, r3, 0.0),), ),
                name='cboom_section1')
            model.rootAssembly.Set(edges=
            model.rootAssembly.instances['cboom-1'].edges.findAt(
                ((0.0, r3, l3),), ),
                name='cboom_section2')
            ##Interaction##
            ##create coupling RP_with_cboomsection
            model.Coupling(name='RP-1_with_cboomsection1', controlPoint=
            Region(referencePoints=(model.rootAssembly.referencePoints[4],)),
                           surface=model.rootAssembly.sets['cboom_section1'],
                           influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC,
                           localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)
            model.Coupling(name='RP-2_with_cboomsection2', controlPoint=
            Region(referencePoints=(model.rootAssembly.referencePoints[5],)),
                           surface=model.rootAssembly.sets['cboom_section2'],
                           influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC,
                           localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)
            ##mesh##
            model.parts[partName].seedPart(deviationFactor=0.1, minSizeFactor=
            0.1, size=meshsize)
            model.parts[partName].generateMesh()
            model.parts[partName].setElementType(elemTypes=(
                ElemType(elemCode=S4R, elemLibrary=EXPLICIT, secondOrderAccuracy=OFF, hourglassControl=STIFFNESS),
                ElemType(elemCode=S3R, elemLibrary=EXPLICIT)), regions=(model.parts[partName].faces,))
            ##set node
            model.rootAssembly.Set(referencePoints=(model.rootAssembly.referencePoints[4],), name='Ref1')
            model.rootAssembly.Set(referencePoints=(model.rootAssembly.referencePoints[5],), name='Ref2')


        print("+" + "-" * 40 + "+")
        print("|{:^40}|".format("Modeling Parameters"))
        print("+" + "-" * 40 + "+")
        print("| {:<12}: {:<25} |".format("Design_space", Design_space))
        print("| {:<12}: {:<25} |".format("length", length))
        print("| {:<12}: {:<25} |".format("radius", radius))
        print("| {:<12}: {:<25} |".format("angle", angle))
        print("+" + "-" * 40 + "+")
        print("+" + "-" * 60 + "+")
        print("|{:^60}|".format("Material Setting Parameters"))
        print("+" + "-" * 60 + "+")
        print("| {0:<15}: {1:<40} |".format("Density", materialdensity))
        print("| {0:<15}: {1:<40} |".format("Material Table", str(materialtable)))
        print("| {0:<15}: {1:<40} |".format("Material Section", str(materialsection)))
        print("+" + "-" * 60 + "+")


    def Static_analysis(Staticnlgeom, Bendingdirection, StaticCPUnum, Staticsaveresult):
        ##Step##
        model.rootAssembly.regenerate()

        if Staticnlgeom == True:
            model.StaticStep(name='Step-bending', previous='Initial', initialInc=0.01, minInc=1e-08, maxInc=0.1,
                             nlgeom=ON)
        else:
            model.StaticStep(name='Step-bending', previous='Initial', initialInc=0.01, minInc=1e-08, maxInc=0.1,
                             nlgeom=OFF)

        model.fieldOutputRequests['F-Output-1'].setValues(variables=(
            'S', 'U', 'RM'))
        model.HistoryOutputRequest(name='H-Output-1', createStepName='Step-bending', variables=('RM2', 'UR2'),
                                   region=model.rootAssembly.sets[
                                       'Ref2'],
                                   sectionPoints=DEFAULT, rebar=EXCLUDE)
        ##Load##
        model.DisplacementBC(name='BC-1', createStepName='Initial',
                             region=model.rootAssembly.sets['Ref1'], u1=SET, u2=SET, u3=SET,
                             ur1=SET, ur2=SET, ur3=SET,
                             amplitude=UNSET, distributionType=UNIFORM, fieldName='',
                             localCsys=None)
        if Bendingdirection == 'X':
            model.DisplacementBC(name='BC-2', createStepName='Step-bending',
                                 region=model.rootAssembly.sets['Ref2'], u1=0.0, u2=UNSET, u3=UNSET,
                                 ur1=-1.0, ur2=0.0, ur3=0.0,
                                 amplitude=UNSET, distributionType=UNIFORM, fieldName='',
                                 localCsys=None)
        if Bendingdirection == 'Y':
            model.DisplacementBC(name='BC-2', createStepName='Step-bending',
                                 region=model.rootAssembly.sets['Ref2'], u1=UNSET, u2=0.0, u3=UNSET,
                                 ur1=0.0, ur2=-1.0, ur3=0.0,
                                 amplitude=UNSET, distributionType=UNIFORM, fieldName='',
                                 localCsys=None)

        ##Job##
        mdb.Job(model='Model-1', name='Static_analysis', numCpus=StaticCPUnum, numDomains=StaticCPUnum)
        mdb.saveAs(pathName='Static_analysis')
        mdb.jobs['Static_analysis'].submit()
        mdb.jobs['Static_analysis'].waitForCompletion()
        mdb.save()
        if Staticsaveresult == True:
            current_directory = os.getcwd()
            odb_filename = "Static_analysis.odb"
            odb_path = os.path.join(current_directory, odb_filename)
            odb = session.openOdb(odb_path)

            session.viewports['Viewport: 1'].setValues(displayedObject=odb)
            odbName = session.viewports[session.currentViewportName].odbDisplay.name

            session.odbData[odbName].setValues(activeFrames=(('Step-bending', (10000,)),))
            session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(
                ('RM', NODAL, ((INVARIANT, 'Magnitude'),)), ('UR', NODAL, ((INVARIANT, 'Magnitude'),)),),
                                        nodeSets=("Ref2",))

            def create_or_open_csv(csv_path):
                if not os.path.exists(csv_path):
                    with open(csv_path, 'w') as file:
                        writer = csv.writer(file)
                        writer.writerow(['UR', 'RM'])
                    print("CSV file created:", csv_path)

            def get_combined_data():
                xy1 = session.xyDataObjects['UR:Magnitude PI: ASSEMBLY N: 2']
                xy2 = session.xyDataObjects['RM:Magnitude PI: ASSEMBLY N: 2']
                xy3 = combine(xy1, xy2)
                data = xy3.data
                return data

            def save_data_to_csv(data, csv_path):
                with open(csv_path, 'wb') as file:
                    writer = csv.writer(file)
                    writer.writerows(data)
                print("Data appended to CSV file:", csv_path)

            data = get_combined_data()
            csv_path = 'UR_RM.csv'
            create_or_open_csv(csv_path)
            save_data_to_csv(data, csv_path)
            del session.xyDataObjects['RM:Magnitude PI: ASSEMBLY N: 2']
            del session.xyDataObjects['UR:Magnitude PI: ASSEMBLY N: 2']
            print("The result has been saved.")

        print("Static analysis has been completed.")

    def NonLinearBuckling_analysis(Bucklingdirection, Bucklingmode, defectcoefficient, saveresult, NonLinearCPUnum):
        model.rootAssembly.regenerate()
        mdb.models.changeKey(fromName='Model-1', toName='Model-LinearBuckling')
        ##Step##
        mdb.models['Model-LinearBuckling'].BuckleStep(name='Step-LinearBuckle', previous='Initial',
                                                      numEigen=20, eigensolver=LANCZOS, minEigen=0.0, blockSize=DEFAULT,
                                                      maxBlocks=DEFAULT)
        mdb.models['Model-LinearBuckling'].fieldOutputRequests['F-Output-1'].setValues(variables=(
            'S', 'U', 'RM'))
        ##Load##
        mdb.models['Model-LinearBuckling'].DisplacementBC(name='BC-1', createStepName='Initial',
                                                          region=mdb.models['Model-LinearBuckling'].rootAssembly.sets['Ref1'], u1=SET, u2=SET,
                                                          u3=SET, ur1=SET,
                                                          ur2=SET, ur3=SET,
                                                          amplitude=UNSET, distributionType=UNIFORM, fieldName='',
                                                          localCsys=None)
        mdb.models['Model-LinearBuckling'].DisplacementBC(name='BC-2', createStepName='Step-LinearBuckle',
                                                          region=mdb.models['Model-LinearBuckling'].rootAssembly.sets['Ref2'], u1=UNSET, u2=SET,
                                                          u3=-1.0, ur1=UNSET,
                                                          ur2=UNSET, ur3=UNSET,
                                                          amplitude=UNSET, distributionType=UNIFORM, fieldName='',
                                                          localCsys=None)
        ##Add keywords
        mdb.models['Model-LinearBuckling'].keywordBlock.synchVersions()
        lines = list(mdb.models['Model-LinearBuckling'].keywordBlock.sieBlocks)  
        target_line = next(i for i, line in enumerate(lines) if "*Restart, write, frequency=0" in line)  
        mdb.models['Model-LinearBuckling'].keywordBlock.insert(target_line , "*node file\n u") 
        ##Job##
        jobName_Linear = 'DCBLinear'
        mdb.Job(model='Model-LinearBuckling', name=jobName_Linear, numCpus=NonLinearCPUnum, numDomains=NonLinearCPUnum)
        mdb.saveAs(pathName=jobName_Linear)
        mdb.jobs[jobName_Linear].submit()
        mdb.jobs[jobName_Linear].waitForCompletion()

        ##copy model
        mdb.Model(name='Model-NonLinearBuckling', objectToCopy=mdb.models['Model-LinearBuckling'])
        del mdb.models['Model-NonLinearBuckling'].steps['Step-LinearBuckle']

        ##NonLinearStep
        mdb.models['Model-NonLinearBuckling'].StaticRiksStep(name='Step-NonLinearBuckle',
                                                             previous='Initial', initialArcInc=0.01, minArcInc=1e-08,
                                                             maxArcInc=0.05,
                                                             nlgeom=ON)
        mdb.models['Model-NonLinearBuckling'].fieldOutputRequests['F-Output-1'].setValues(variables=('S', 'U', 'RM'))
        mdb.models['Model-NonLinearBuckling'].FieldOutputRequest(name='F-Output-2',
                                                                 createStepName='Step-NonLinearBuckle',
                                                                 variables=('S', 'U', 'UR', 'RM'),
                                                                 region=mdb.models[
                                                                     'Model-NonLinearBuckling'].rootAssembly.sets[
                                                                     'Ref2'],
                                                                 sectionPoints=DEFAULT, rebar=EXCLUDE)
        mdb.models['Model-NonLinearBuckling'].HistoryOutputRequest(name='H-Output-2',
                                                                   createStepName='Step-NonLinearBuckle',
                                                                   variables=('UR', 'RM'),
                                                                   region=mdb.models[
                                                                       'Model-NonLinearBuckling'].rootAssembly.sets[
                                                                       'Ref2'], sectionPoints=DEFAULT, rebar=EXCLUDE)
        ##NonLinearLoad
        mdb.models['Model-NonLinearBuckling'].DisplacementBC(name='BC-1', createStepName='Initial',
                                                             region=mdb.models[
                                                                       'Model-NonLinearBuckling'].rootAssembly.sets['Ref1'], u1=SET, u2=SET,
                                                             u3=SET, ur1=SET,
                                                             ur2=SET, ur3=SET,
                                                             amplitude=UNSET, distributionType=UNIFORM, fieldName='',
                                                             localCsys=None)
        if Bucklingdirection == 'X':
            mdb.models['Model-NonLinearBuckling'].DisplacementBC(name='BC-2', createStepName='Step-NonLinearBuckle',
                                                                 region=mdb.models[
                                                                       'Model-NonLinearBuckling'].rootAssembly.sets['Ref2'], u1=0.0,
                                                                 u2=UNSET, u3=UNSET,
                                                                 ur1=-1.0, ur2=0.0, ur3=0.0,
                                                                 amplitude=UNSET, distributionType=UNIFORM,
                                                                 fieldName='',
                                                                 localCsys=None)
        if Bucklingdirection == 'Y':
            mdb.models['Model-NonLinearBuckling'].DisplacementBC(name='BC-2', createStepName='Step-NonLinearBuckle',
                                                                 region=mdb.models[
                                                                       'Model-NonLinearBuckling'].rootAssembly.sets['Ref2'], u1=UNSET,
                                                                 u2=0.0, u3=UNSET,
                                                                 ur1=0.0, ur2=-1.0, ur3=0.0,
                                                                 amplitude=UNSET, distributionType=UNIFORM,
                                                                 fieldName='',
                                                                 localCsys=None)

        ##Add keywords

        mdb.models['Model-NonLinearBuckling'].keywordBlock.setValues(edited=0)
        mdb.models['Model-NonLinearBuckling'].keywordBlock.synchVersions(
        storeNodesAndElements=False)
        lines = list(mdb.models['Model-NonLinearBuckling'].keywordBlock.sieBlocks)
        target_line = next(i for i, line in enumerate(lines) if "** STEP: Step-NonLinearBuckle" in line)
        new_line = "*imperfection,file=DCBLinear,Step=1\n{},{}\n".format(Bucklingmode, defectcoefficient)
        mdb.models['Model-NonLinearBuckling'].keywordBlock.insert(target_line-1, new_line)

        ##Job##
        jobName_NonLinear = 'DCBNonLinear'
        mdb.Job(model='Model-NonLinearBuckling', name=jobName_NonLinear, numCpus=NonLinearCPUnum,
                numDomains=NonLinearCPUnum)
        mdb.saveAs(pathName=jobName_NonLinear)
        mdb.jobs[jobName_NonLinear].submit()
        mdb.jobs[jobName_NonLinear].waitForCompletion()
        mdb.save()
        if saveresult == True:
            current_directory = os.getcwd()
            odb_filename = "DCBNonLinear.odb"
            odb_path = os.path.join(current_directory, odb_filename)
            odb = session.openOdb(odb_path)

            session.viewports['Viewport: 1'].setValues(displayedObject=odb)
            odbName = session.viewports[session.currentViewportName].odbDisplay.name

            session.odbData[odbName].setValues(activeFrames=(('Step-NonLinearBuckle', (10000,)),))
            session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(
                ('RM', NODAL, ((INVARIANT, 'Magnitude'),)), ('UR', NODAL, ((INVARIANT, 'Magnitude'),)),),
                                        nodeSets=("Ref2",))

            def create_or_open_csv(csv_path):
                if not os.path.exists(csv_path):
                    with open(csv_path, 'w') as file:
                        writer = csv.writer(file)
                        writer.writerow(['UR', 'RM'])
                    print("CSV file created:", csv_path)

            def get_combined_data():
                xy1 = session.xyDataObjects['UR:Magnitude PI: ASSEMBLY N: 2']
                xy2 = session.xyDataObjects['RM:Magnitude PI: ASSEMBLY N: 2']
                xy3 = combine(xy1, xy2)
                data = xy3.data
                return data

            def save_data_to_csv(data, csv_path):
                with open(csv_path, 'w') as file:
                    writer = csv.writer(file)
                    writer.writerows(data)
                print("Data appended to CSV file:", csv_path)

            data = get_combined_data()
            csv_path = 'UR_RM.csv'
            create_or_open_csv(csv_path)
            save_data_to_csv(data, csv_path)
            del session.xyDataObjects['RM:Magnitude PI: ASSEMBLY N: 2']
            del session.xyDataObjects['UR:Magnitude PI: ASSEMBLY N: 2']
            print("The result has been saved.")

        print("Riks method-based nonlinear buckling analysis has been completed.")

    def Frequency_analysis(Eigenreq, RefBC, Inertiamass, freqCPUnum):
        ##Step##
        model.rootAssembly.regenerate()
        model.FrequencyStep(name='Step-freq', previous='Initial',
                            numEigen=Eigenreq)
        model.fieldOutputRequests['F-Output-1'].setValues(
            variables=('S', 'U'))
        ##Load##
        if RefBC == True:
            model.DisplacementBC(name='BC-1', createStepName='Initial',
                                 region=model.rootAssembly.sets['Ref1'], u1=SET,
                                 u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET,
                                 amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)

        model.rootAssembly.engineeringFeatures.PointMassInertia(
            name='Inertia-1', region=model.rootAssembly.sets['Ref2'], mass=Inertiamass, alpha=0.0, composite=0.0)

        ##Job##
        FreqjobName = 'DCBfreq'
        mdb.Job(model='Model-1', name=FreqjobName, numCpus=freqCPUnum, numDomains=freqCPUnum)
        mdb.saveAs(pathName=FreqjobName)
        mdb.jobs[FreqjobName].submit()
        mdb.jobs[FreqjobName].waitForCompletion()
        mdb.save()

        print("Frequency analysis has been completed.")


    Modeling(model_lib, Design_space, length, radius, angle,materialdensity, materialtable, materialsection,meshsize)


    if   anslysislib == 'No Analysis':
         print('Do not conduct analysis')

    elif anslysislib == 'Static Analysis':
         print('Static Analysis')
         Static_analysis(Staticnlgeom, Bendingdirection, StaticCPUnum, Staticsaveresult)

    elif anslysislib == 'NonLinearBuckling Analysis':
         print('NonLinearBuckling Analysis')
         NonLinearBuckling_analysis(Bucklingdirection, Bucklingmode, defectcoefficient, saveresult, NonLinearCPUnum)

    elif anslysislib == 'Frequency Analysis':
         print('Frequency Analysis')
         Frequency_analysis(Eigenreq, RefBC, Inertiamass, freqCPUnum)

