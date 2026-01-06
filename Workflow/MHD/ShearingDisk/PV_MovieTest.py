# trace generated using paraview version 5.10.1
import paraview
import glob
import sys
import os
import PV_VarDict

paraview.compatibility.major = 5
paraview.compatibility.minor = 10

from PV_VarDict import SD_Params

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'AMReX/BoxLib Grid Reader'

PlotFileDir = sys.argv[1]
PlotFileList = sorted(glob.glob( PlotFileDir + '/' + 'ShearingDisk2D.plt.*' ))
PlotVar = sys.argv[2]
PlotParams = SD_Params[PlotVar]

shearingDisk2Dplt0 = AMReXBoxLibGridReader(registrationName='ShearingDisk2D.plt.0*', FileNames = PlotFileList )
shearingDisk2Dplt0.EnableCaching = 0
shearingDisk2Dplt0.Level = 1
shearingDisk2Dplt0.PointArrayStatus = []
shearingDisk2Dplt0.CellArrayStatus = []

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# Properties modified on shearingDisk2Dplt0
shearingDisk2Dplt0.CellArrayStatus = ['CM_B1', 'CM_B2', 'CM_B3', 'PM_V1', 'PM_V3', 'CM_D']

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
shearingDisk2Dplt0Display = Show(shearingDisk2Dplt0, renderView1, 'AMRRepresentation')

# trace defaults for the display properties.
shearingDisk2Dplt0Display.Selection = None
shearingDisk2Dplt0Display.Representation = 'Outline'
shearingDisk2Dplt0Display.ColorArrayName = [None, '']
shearingDisk2Dplt0Display.LookupTable = None
shearingDisk2Dplt0Display.MapScalars = 1
shearingDisk2Dplt0Display.MultiComponentsMapping = 0
shearingDisk2Dplt0Display.InterpolateScalarsBeforeMapping = 1
shearingDisk2Dplt0Display.Opacity = 1.0
shearingDisk2Dplt0Display.PointSize = 2.0
shearingDisk2Dplt0Display.LineWidth = 1.0
shearingDisk2Dplt0Display.RenderLinesAsTubes = 0
shearingDisk2Dplt0Display.RenderPointsAsSpheres = 0
shearingDisk2Dplt0Display.Interpolation = 'Gouraud'
shearingDisk2Dplt0Display.Specular = 0.0
shearingDisk2Dplt0Display.SpecularColor = [1.0, 1.0, 1.0]
shearingDisk2Dplt0Display.SpecularPower = 100.0
shearingDisk2Dplt0Display.Luminosity = 0.0
shearingDisk2Dplt0Display.Ambient = 0.0
shearingDisk2Dplt0Display.Diffuse = 1.0
shearingDisk2Dplt0Display.Roughness = 0.3
shearingDisk2Dplt0Display.Metallic = 0.0
shearingDisk2Dplt0Display.EdgeTint = [1.0, 1.0, 1.0]
shearingDisk2Dplt0Display.Anisotropy = 0.0
shearingDisk2Dplt0Display.AnisotropyRotation = 0.0
shearingDisk2Dplt0Display.BaseIOR = 1.5
shearingDisk2Dplt0Display.CoatStrength = 0.0
shearingDisk2Dplt0Display.CoatIOR = 2.0
shearingDisk2Dplt0Display.CoatRoughness = 0.0
shearingDisk2Dplt0Display.CoatColor = [1.0, 1.0, 1.0]
shearingDisk2Dplt0Display.SelectTCoordArray = 'None'
shearingDisk2Dplt0Display.SelectNormalArray = 'None'
shearingDisk2Dplt0Display.SelectTangentArray = 'None'
shearingDisk2Dplt0Display.Texture = None
shearingDisk2Dplt0Display.RepeatTextures = 1
shearingDisk2Dplt0Display.InterpolateTextures = 0
shearingDisk2Dplt0Display.SeamlessU = 0
shearingDisk2Dplt0Display.SeamlessV = 0
shearingDisk2Dplt0Display.UseMipmapTextures = 0
shearingDisk2Dplt0Display.ShowTexturesOnBackface = 1
shearingDisk2Dplt0Display.BaseColorTexture = None
shearingDisk2Dplt0Display.NormalTexture = None
shearingDisk2Dplt0Display.NormalScale = 1.0
shearingDisk2Dplt0Display.CoatNormalTexture = None
shearingDisk2Dplt0Display.CoatNormalScale = 1.0
shearingDisk2Dplt0Display.MaterialTexture = None
shearingDisk2Dplt0Display.OcclusionStrength = 1.0
shearingDisk2Dplt0Display.AnisotropyTexture = None
shearingDisk2Dplt0Display.EmissiveTexture = None
shearingDisk2Dplt0Display.EmissiveFactor = [1.0, 1.0, 1.0]
shearingDisk2Dplt0Display.FlipTextures = 0
shearingDisk2Dplt0Display.BackfaceRepresentation = 'Follow Frontface'
shearingDisk2Dplt0Display.BackfaceAmbientColor = [1.0, 1.0, 1.0]
shearingDisk2Dplt0Display.BackfaceOpacity = 1.0
shearingDisk2Dplt0Display.Position = [0.0, 0.0, 0.0]
shearingDisk2Dplt0Display.Scale = [1.0, 1.0, 1.0]
shearingDisk2Dplt0Display.Orientation = [0.0, 0.0, 0.0]
shearingDisk2Dplt0Display.Origin = [0.0, 0.0, 0.0]
shearingDisk2Dplt0Display.CoordinateShiftScaleMethod = 'Always Auto Shift Scale'
shearingDisk2Dplt0Display.Pickable = 1
shearingDisk2Dplt0Display.Triangulate = 0
shearingDisk2Dplt0Display.UseShaderReplacements = 0
shearingDisk2Dplt0Display.ShaderReplacements = ''
shearingDisk2Dplt0Display.NonlinearSubdivisionLevel = 1
shearingDisk2Dplt0Display.UseDataPartitions = 0
shearingDisk2Dplt0Display.OSPRayUseScaleArray = 'All Approximate'
shearingDisk2Dplt0Display.OSPRayScaleArray = ''
shearingDisk2Dplt0Display.OSPRayScaleFunction = 'PiecewiseFunction'
shearingDisk2Dplt0Display.OSPRayMaterial = 'None'
shearingDisk2Dplt0Display.BlockSelectors = ['/']
shearingDisk2Dplt0Display.BlockColors = []
shearingDisk2Dplt0Display.BlockOpacities = []
shearingDisk2Dplt0Display.Orient = 0
shearingDisk2Dplt0Display.OrientationMode = 'Direction'
shearingDisk2Dplt0Display.SelectOrientationVectors = 'None'
shearingDisk2Dplt0Display.Scaling = 0
shearingDisk2Dplt0Display.ScaleMode = 'No Data Scaling Off'
shearingDisk2Dplt0Display.ScaleFactor = 0.1
shearingDisk2Dplt0Display.SelectScaleArray = 'None'
shearingDisk2Dplt0Display.GlyphType = 'Arrow'
shearingDisk2Dplt0Display.UseGlyphTable = 0
shearingDisk2Dplt0Display.GlyphTableIndexArray = 'None'
shearingDisk2Dplt0Display.UseCompositeGlyphTable = 0
shearingDisk2Dplt0Display.UseGlyphCullingAndLOD = 0
shearingDisk2Dplt0Display.LODValues = []
shearingDisk2Dplt0Display.ColorByLODIndex = 0
shearingDisk2Dplt0Display.GaussianRadius = 0.005
shearingDisk2Dplt0Display.ShaderPreset = 'Sphere'
shearingDisk2Dplt0Display.CustomTriangleScale = 3
shearingDisk2Dplt0Display.CustomShader = """ // This custom shader code define a gaussian blur
 // Please take a look into vtkSMPointGaussianRepresentation.cxx
 // for other custom shader examples
 //VTK::Color::Impl
   float dist2 = dot(offsetVCVSOutput.xy,offsetVCVSOutput.xy);
   float gaussian = exp(-0.5*dist2);
   opacity = opacity*gaussian;
"""
shearingDisk2Dplt0Display.Emissive = 0
shearingDisk2Dplt0Display.ScaleByArray = 0
shearingDisk2Dplt0Display.SetScaleArray = [None, '']
shearingDisk2Dplt0Display.ScaleArrayComponent = 0
shearingDisk2Dplt0Display.UseScaleFunction = 1
shearingDisk2Dplt0Display.ScaleTransferFunction = 'PiecewiseFunction'
shearingDisk2Dplt0Display.OpacityByArray = 0
shearingDisk2Dplt0Display.OpacityArray = [None, '']
shearingDisk2Dplt0Display.OpacityArrayComponent = 0
shearingDisk2Dplt0Display.OpacityTransferFunction = 'PiecewiseFunction'
shearingDisk2Dplt0Display.DataAxesGrid = 'GridAxesRepresentation'
shearingDisk2Dplt0Display.SelectionCellLabelBold = 0
shearingDisk2Dplt0Display.SelectionCellLabelColor = [0.0, 1.0, 0.0]
shearingDisk2Dplt0Display.SelectionCellLabelFontFamily = 'Arial'
shearingDisk2Dplt0Display.SelectionCellLabelFontFile = ''
shearingDisk2Dplt0Display.SelectionCellLabelFontSize = 24
shearingDisk2Dplt0Display.SelectionCellLabelItalic = 0
shearingDisk2Dplt0Display.SelectionCellLabelJustification = 'Left'
shearingDisk2Dplt0Display.SelectionCellLabelOpacity = 1.0
shearingDisk2Dplt0Display.SelectionCellLabelShadow = 0
shearingDisk2Dplt0Display.SelectionPointLabelBold = 1
shearingDisk2Dplt0Display.SelectionPointLabelColor = [1.0, 1.0, 0.0]
shearingDisk2Dplt0Display.SelectionPointLabelFontFamily = 'Arial'
shearingDisk2Dplt0Display.SelectionPointLabelFontFile = ''
shearingDisk2Dplt0Display.SelectionPointLabelFontSize = 18
shearingDisk2Dplt0Display.SelectionPointLabelItalic = 0
shearingDisk2Dplt0Display.SelectionPointLabelJustification = 'Left'
shearingDisk2Dplt0Display.SelectionPointLabelOpacity = 1.0
shearingDisk2Dplt0Display.SelectionPointLabelShadow = 0
shearingDisk2Dplt0Display.PolarAxes = 'PolarAxesRepresentation'
shearingDisk2Dplt0Display.ScalarOpacityUnitDistance = 0.06564197879454707
shearingDisk2Dplt0Display.ScalarOpacityFunction = None
shearingDisk2Dplt0Display.VolumeRenderingMode = 'Smart'
shearingDisk2Dplt0Display.ResamplingMode = 'Over Data Bounds'
shearingDisk2Dplt0Display.StreamingRequestSize = 10
shearingDisk2Dplt0Display.NumberOfSamples = [64, 128, 64]
shearingDisk2Dplt0Display.Shade = 0

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
shearingDisk2Dplt0Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
shearingDisk2Dplt0Display.OSPRayScaleFunction.UseLogScale = 0

# init the 'Arrow' selected for 'GlyphType'
shearingDisk2Dplt0Display.GlyphType.TipResolution = 6
shearingDisk2Dplt0Display.GlyphType.TipRadius = 0.1
shearingDisk2Dplt0Display.GlyphType.TipLength = 0.35
shearingDisk2Dplt0Display.GlyphType.ShaftResolution = 6
shearingDisk2Dplt0Display.GlyphType.ShaftRadius = 0.03
shearingDisk2Dplt0Display.GlyphType.Invert = 0

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
shearingDisk2Dplt0Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
shearingDisk2Dplt0Display.ScaleTransferFunction.UseLogScale = 0

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
shearingDisk2Dplt0Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
shearingDisk2Dplt0Display.OpacityTransferFunction.UseLogScale = 0

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
shearingDisk2Dplt0Display.DataAxesGrid.XTitle = 'X Axis'
shearingDisk2Dplt0Display.DataAxesGrid.YTitle = 'Y Axis'
shearingDisk2Dplt0Display.DataAxesGrid.ZTitle = 'Z Axis'
shearingDisk2Dplt0Display.DataAxesGrid.XTitleFontFamily = 'Arial'
shearingDisk2Dplt0Display.DataAxesGrid.XTitleFontFile = ''
shearingDisk2Dplt0Display.DataAxesGrid.XTitleBold = 0
shearingDisk2Dplt0Display.DataAxesGrid.XTitleItalic = 0
shearingDisk2Dplt0Display.DataAxesGrid.XTitleFontSize = 12
shearingDisk2Dplt0Display.DataAxesGrid.XTitleShadow = 0
shearingDisk2Dplt0Display.DataAxesGrid.XTitleOpacity = 1.0
shearingDisk2Dplt0Display.DataAxesGrid.YTitleFontFamily = 'Arial'
shearingDisk2Dplt0Display.DataAxesGrid.YTitleFontFile = ''
shearingDisk2Dplt0Display.DataAxesGrid.YTitleBold = 0
shearingDisk2Dplt0Display.DataAxesGrid.YTitleItalic = 0
shearingDisk2Dplt0Display.DataAxesGrid.YTitleFontSize = 12
shearingDisk2Dplt0Display.DataAxesGrid.YTitleShadow = 0
shearingDisk2Dplt0Display.DataAxesGrid.YTitleOpacity = 1.0
shearingDisk2Dplt0Display.DataAxesGrid.ZTitleFontFamily = 'Arial'
shearingDisk2Dplt0Display.DataAxesGrid.ZTitleFontFile = ''
shearingDisk2Dplt0Display.DataAxesGrid.ZTitleBold = 0
shearingDisk2Dplt0Display.DataAxesGrid.ZTitleItalic = 0
shearingDisk2Dplt0Display.DataAxesGrid.ZTitleFontSize = 12
shearingDisk2Dplt0Display.DataAxesGrid.ZTitleShadow = 0
shearingDisk2Dplt0Display.DataAxesGrid.ZTitleOpacity = 1.0
shearingDisk2Dplt0Display.DataAxesGrid.FacesToRender = 63
shearingDisk2Dplt0Display.DataAxesGrid.CullBackface = 0
shearingDisk2Dplt0Display.DataAxesGrid.CullFrontface = 1
shearingDisk2Dplt0Display.DataAxesGrid.ShowGrid = 0
shearingDisk2Dplt0Display.DataAxesGrid.ShowEdges = 1
shearingDisk2Dplt0Display.DataAxesGrid.ShowTicks = 1
shearingDisk2Dplt0Display.DataAxesGrid.LabelUniqueEdgesOnly = 1
shearingDisk2Dplt0Display.DataAxesGrid.AxesToLabel = 63
shearingDisk2Dplt0Display.DataAxesGrid.XLabelFontFamily = 'Arial'
shearingDisk2Dplt0Display.DataAxesGrid.XLabelFontFile = ''
shearingDisk2Dplt0Display.DataAxesGrid.XLabelBold = 0
shearingDisk2Dplt0Display.DataAxesGrid.XLabelItalic = 0
shearingDisk2Dplt0Display.DataAxesGrid.XLabelFontSize = 12
shearingDisk2Dplt0Display.DataAxesGrid.XLabelShadow = 0
shearingDisk2Dplt0Display.DataAxesGrid.XLabelOpacity = 1.0
shearingDisk2Dplt0Display.DataAxesGrid.YLabelFontFamily = 'Arial'
shearingDisk2Dplt0Display.DataAxesGrid.YLabelFontFile = ''
shearingDisk2Dplt0Display.DataAxesGrid.YLabelBold = 0
shearingDisk2Dplt0Display.DataAxesGrid.YLabelItalic = 0
shearingDisk2Dplt0Display.DataAxesGrid.YLabelFontSize = 12
shearingDisk2Dplt0Display.DataAxesGrid.YLabelShadow = 0
shearingDisk2Dplt0Display.DataAxesGrid.YLabelOpacity = 1.0
shearingDisk2Dplt0Display.DataAxesGrid.ZLabelFontFamily = 'Arial'
shearingDisk2Dplt0Display.DataAxesGrid.ZLabelFontFile = ''
shearingDisk2Dplt0Display.DataAxesGrid.ZLabelBold = 0
shearingDisk2Dplt0Display.DataAxesGrid.ZLabelItalic = 0
shearingDisk2Dplt0Display.DataAxesGrid.ZLabelFontSize = 12
shearingDisk2Dplt0Display.DataAxesGrid.ZLabelShadow = 0
shearingDisk2Dplt0Display.DataAxesGrid.ZLabelOpacity = 1.0
shearingDisk2Dplt0Display.DataAxesGrid.XAxisNotation = 'Mixed'
shearingDisk2Dplt0Display.DataAxesGrid.XAxisPrecision = 2
shearingDisk2Dplt0Display.DataAxesGrid.XAxisUseCustomLabels = 0
shearingDisk2Dplt0Display.DataAxesGrid.XAxisLabels = []
shearingDisk2Dplt0Display.DataAxesGrid.YAxisNotation = 'Mixed'
shearingDisk2Dplt0Display.DataAxesGrid.YAxisPrecision = 2
shearingDisk2Dplt0Display.DataAxesGrid.YAxisUseCustomLabels = 0
shearingDisk2Dplt0Display.DataAxesGrid.YAxisLabels = []
shearingDisk2Dplt0Display.DataAxesGrid.ZAxisNotation = 'Mixed'
shearingDisk2Dplt0Display.DataAxesGrid.ZAxisPrecision = 2
shearingDisk2Dplt0Display.DataAxesGrid.ZAxisUseCustomLabels = 0
shearingDisk2Dplt0Display.DataAxesGrid.ZAxisLabels = []
shearingDisk2Dplt0Display.DataAxesGrid.UseCustomBounds = 0
shearingDisk2Dplt0Display.DataAxesGrid.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
shearingDisk2Dplt0Display.PolarAxes.Visibility = 0
shearingDisk2Dplt0Display.PolarAxes.Translation = [0.0, 0.0, 0.0]
shearingDisk2Dplt0Display.PolarAxes.Scale = [1.0, 1.0, 1.0]
shearingDisk2Dplt0Display.PolarAxes.Orientation = [0.0, 0.0, 0.0]
shearingDisk2Dplt0Display.PolarAxes.EnableCustomBounds = [0, 0, 0]
shearingDisk2Dplt0Display.PolarAxes.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
shearingDisk2Dplt0Display.PolarAxes.EnableCustomRange = 0
shearingDisk2Dplt0Display.PolarAxes.CustomRange = [0.0, 1.0]
shearingDisk2Dplt0Display.PolarAxes.PolarAxisVisibility = 1
shearingDisk2Dplt0Display.PolarAxes.RadialAxesVisibility = 1
shearingDisk2Dplt0Display.PolarAxes.DrawRadialGridlines = 1
shearingDisk2Dplt0Display.PolarAxes.PolarArcsVisibility = 1
shearingDisk2Dplt0Display.PolarAxes.DrawPolarArcsGridlines = 1
shearingDisk2Dplt0Display.PolarAxes.NumberOfRadialAxes = 0
shearingDisk2Dplt0Display.PolarAxes.AutoSubdividePolarAxis = 1
shearingDisk2Dplt0Display.PolarAxes.NumberOfPolarAxis = 0
shearingDisk2Dplt0Display.PolarAxes.MinimumRadius = 0.0
shearingDisk2Dplt0Display.PolarAxes.MinimumAngle = 0.0
shearingDisk2Dplt0Display.PolarAxes.MaximumAngle = 90.0
shearingDisk2Dplt0Display.PolarAxes.RadialAxesOriginToPolarAxis = 1
shearingDisk2Dplt0Display.PolarAxes.Ratio = 1.0
shearingDisk2Dplt0Display.PolarAxes.PolarAxisColor = [1.0, 1.0, 1.0]
shearingDisk2Dplt0Display.PolarAxes.PolarArcsColor = [1.0, 1.0, 1.0]
shearingDisk2Dplt0Display.PolarAxes.LastRadialAxisColor = [1.0, 1.0, 1.0]
shearingDisk2Dplt0Display.PolarAxes.SecondaryPolarArcsColor = [1.0, 1.0, 1.0]
shearingDisk2Dplt0Display.PolarAxes.SecondaryRadialAxesColor = [1.0, 1.0, 1.0]
shearingDisk2Dplt0Display.PolarAxes.PolarAxisTitleVisibility = 1
shearingDisk2Dplt0Display.PolarAxes.PolarAxisTitle = 'Radial Distance'
shearingDisk2Dplt0Display.PolarAxes.PolarAxisTitleLocation = 'Bottom'
shearingDisk2Dplt0Display.PolarAxes.PolarLabelVisibility = 1
shearingDisk2Dplt0Display.PolarAxes.PolarLabelFormat = '%-#6.3g'
shearingDisk2Dplt0Display.PolarAxes.PolarLabelExponentLocation = 'Labels'
shearingDisk2Dplt0Display.PolarAxes.RadialLabelVisibility = 1
shearingDisk2Dplt0Display.PolarAxes.RadialLabelFormat = '%-#3.1f'
shearingDisk2Dplt0Display.PolarAxes.RadialLabelLocation = 'Bottom'
shearingDisk2Dplt0Display.PolarAxes.RadialUnitsVisibility = 1
shearingDisk2Dplt0Display.PolarAxes.ScreenSize = 10.0
shearingDisk2Dplt0Display.PolarAxes.PolarAxisTitleOpacity = 1.0
shearingDisk2Dplt0Display.PolarAxes.PolarAxisTitleFontFamily = 'Arial'
shearingDisk2Dplt0Display.PolarAxes.PolarAxisTitleFontFile = ''
shearingDisk2Dplt0Display.PolarAxes.PolarAxisTitleBold = 0
shearingDisk2Dplt0Display.PolarAxes.PolarAxisTitleItalic = 0
shearingDisk2Dplt0Display.PolarAxes.PolarAxisTitleShadow = 0
shearingDisk2Dplt0Display.PolarAxes.PolarAxisTitleFontSize = 12
shearingDisk2Dplt0Display.PolarAxes.PolarAxisLabelOpacity = 1.0
shearingDisk2Dplt0Display.PolarAxes.PolarAxisLabelFontFamily = 'Arial'
shearingDisk2Dplt0Display.PolarAxes.PolarAxisLabelFontFile = ''
shearingDisk2Dplt0Display.PolarAxes.PolarAxisLabelBold = 0
shearingDisk2Dplt0Display.PolarAxes.PolarAxisLabelItalic = 0
shearingDisk2Dplt0Display.PolarAxes.PolarAxisLabelShadow = 0
shearingDisk2Dplt0Display.PolarAxes.PolarAxisLabelFontSize = 12
shearingDisk2Dplt0Display.PolarAxes.LastRadialAxisTextOpacity = 1.0
shearingDisk2Dplt0Display.PolarAxes.LastRadialAxisTextFontFamily = 'Arial'
shearingDisk2Dplt0Display.PolarAxes.LastRadialAxisTextFontFile = ''
shearingDisk2Dplt0Display.PolarAxes.LastRadialAxisTextBold = 0
shearingDisk2Dplt0Display.PolarAxes.LastRadialAxisTextItalic = 0
shearingDisk2Dplt0Display.PolarAxes.LastRadialAxisTextShadow = 0
shearingDisk2Dplt0Display.PolarAxes.LastRadialAxisTextFontSize = 12
shearingDisk2Dplt0Display.PolarAxes.SecondaryRadialAxesTextOpacity = 1.0
shearingDisk2Dplt0Display.PolarAxes.SecondaryRadialAxesTextFontFamily = 'Arial'
shearingDisk2Dplt0Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''
shearingDisk2Dplt0Display.PolarAxes.SecondaryRadialAxesTextBold = 0
shearingDisk2Dplt0Display.PolarAxes.SecondaryRadialAxesTextItalic = 0
shearingDisk2Dplt0Display.PolarAxes.SecondaryRadialAxesTextShadow = 0
shearingDisk2Dplt0Display.PolarAxes.SecondaryRadialAxesTextFontSize = 12
shearingDisk2Dplt0Display.PolarAxes.EnableDistanceLOD = 1
shearingDisk2Dplt0Display.PolarAxes.DistanceLODThreshold = 0.7
shearingDisk2Dplt0Display.PolarAxes.EnableViewAngleLOD = 1
shearingDisk2Dplt0Display.PolarAxes.ViewAngleLODThreshold = 0.7
shearingDisk2Dplt0Display.PolarAxes.SmallestVisiblePolarAngle = 0.5
shearingDisk2Dplt0Display.PolarAxes.PolarTicksVisibility = 1
shearingDisk2Dplt0Display.PolarAxes.ArcTicksOriginToPolarAxis = 1
shearingDisk2Dplt0Display.PolarAxes.TickLocation = 'Both'
shearingDisk2Dplt0Display.PolarAxes.AxisTickVisibility = 1
shearingDisk2Dplt0Display.PolarAxes.AxisMinorTickVisibility = 0
shearingDisk2Dplt0Display.PolarAxes.ArcTickVisibility = 1
shearingDisk2Dplt0Display.PolarAxes.ArcMinorTickVisibility = 0
shearingDisk2Dplt0Display.PolarAxes.DeltaAngleMajor = 10.0
shearingDisk2Dplt0Display.PolarAxes.DeltaAngleMinor = 5.0
shearingDisk2Dplt0Display.PolarAxes.PolarAxisMajorTickSize = 0.0
shearingDisk2Dplt0Display.PolarAxes.PolarAxisTickRatioSize = 0.3
shearingDisk2Dplt0Display.PolarAxes.PolarAxisMajorTickThickness = 1.0
shearingDisk2Dplt0Display.PolarAxes.PolarAxisTickRatioThickness = 0.5
shearingDisk2Dplt0Display.PolarAxes.LastRadialAxisMajorTickSize = 0.0
shearingDisk2Dplt0Display.PolarAxes.LastRadialAxisTickRatioSize = 0.3
shearingDisk2Dplt0Display.PolarAxes.LastRadialAxisMajorTickThickness = 1.0
shearingDisk2Dplt0Display.PolarAxes.LastRadialAxisTickRatioThickness = 0.5
shearingDisk2Dplt0Display.PolarAxes.ArcMajorTickSize = 0.0
shearingDisk2Dplt0Display.PolarAxes.ArcTickRatioSize = 0.3
shearingDisk2Dplt0Display.PolarAxes.ArcMajorTickThickness = 1.0
shearingDisk2Dplt0Display.PolarAxes.ArcTickRatioThickness = 0.5
shearingDisk2Dplt0Display.PolarAxes.Use2DMode = 0
shearingDisk2Dplt0Display.PolarAxes.UseLogAxis = 0

# reset view to fit data
renderView1.ResetCamera(False)

#changing interaction mode based on data extents
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [15.5, 0.0, 10000.0]
renderView1.CameraFocalPoint = [15.5, 0.0, 0.0]

# get the material library
materialLibrary1 = GetMaterialLibrary()

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(shearingDisk2Dplt0Display, ('CELLS', PlotVar))

VarDLUT = GetColorTransferFunction( PlotVar )
VarDLUT.ApplyPreset( PlotParams[0] )
VarDLUT.RescaleTransferFunction( PlotParams[1], PlotParams[2] )

# show color bar/color legend
shearingDisk2Dplt0Display.SetScalarBarVisibility(renderView1, True)

# change representation type
shearingDisk2Dplt0Display.SetRepresentationType('Surface')

# create a new 'Annotate Time Filter'
annotateTimeFilter1 = AnnotateTimeFilter(registrationName='AnnotateTimeFilter1', Input=shearingDisk2Dplt0)
annotateTimeFilter1.Format = 'Time: {time:f}'
annotateTimeFilter1.Shift = 0.0
annotateTimeFilter1.Scale = 1.0

# set active source
SetActiveSource(annotateTimeFilter1)

# show data in view
annotateTimeFilter1Display = Show(annotateTimeFilter1, renderView1, 'TextSourceRepresentation')

# trace defaults for the display properties.
annotateTimeFilter1Display.TextPropMode = '2D Text Widget'
annotateTimeFilter1Display.Interactivity = 1
annotateTimeFilter1Display.WindowLocation = 'Upper Left Corner'
annotateTimeFilter1Display.Position = [0.05, 0.05]
annotateTimeFilter1Display.Opacity = 1.0
annotateTimeFilter1Display.FontFamily = 'Arial'
annotateTimeFilter1Display.FontFile = ''
annotateTimeFilter1Display.Bold = 0
annotateTimeFilter1Display.Italic = 0
annotateTimeFilter1Display.Shadow = 0
annotateTimeFilter1Display.FontSize = 18
annotateTimeFilter1Display.Justification = 'Center'
annotateTimeFilter1Display.VerticalJustification = 'Center'
annotateTimeFilter1Display.ShowBorder = 'Only on hover'
annotateTimeFilter1Display.BackgroundColor = [1.0, 1.0, 1.0, 0.2]
annotateTimeFilter1Display.BorderThickness = 0.0
annotateTimeFilter1Display.CornerRadius = 0.0
annotateTimeFilter1Display.Padding = 1
annotateTimeFilter1Display.BasePosition = [0.0, 0.0, 0.0]
annotateTimeFilter1Display.TopPosition = [0.0, 1.0, 0.0]
annotateTimeFilter1Display.FlagSize = 1.0
annotateTimeFilter1Display.BillboardPosition = [0.0, 0.0, 0.0]
annotateTimeFilter1Display.DisplayOffset = [0, 0]

# show data in view
annotateTimeFilter1Display = Show(annotateTimeFilter1, renderView1, 'TextSourceRepresentation')

# get layout
layout1 = GetLayout()

# layout/tab size in pixels
layout1.SetSize(716, 391)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [15.5, 0.0, 10000.0]
renderView1.CameraFocalPoint = [15.5, 0.0, 0.0]
renderView1.CameraParallelScale = 0.7071067811865476

# save animation
SaveAnimation( PlotFileDir + '/' + PlotVar + '_mov.png', renderView1, ImageResolution=[1432, 782],
    FontScaling='Scale fonts proportionally',
    OverrideColorPalette='',
    StereoMode='No change',
    TransparentBackground=0,
    FrameRate=1,
    FrameWindow=[0, 121], 
    # PNG options
    CompressionLevel='0',
    MetaData=['Application', 'ParaView'],
    SuffixFormat='.%04d')

#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(716, 391)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [15.5, 0.0, 10000.0]
renderView1.CameraFocalPoint = [15.5, 0.0, 0.0]
renderView1.CameraParallelScale = 0.7071067811865476

#--------------------------------------------
# uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).

print("Making movie.")
MovCmd = "ffmpeg -f image2 -r 10 -i " + PlotFileDir + "/" + PlotVar + "_mov.%04d.png -vcodec mpeg4 -y " + PlotFileDir + "/" + PlotVar + ".mp4"

os.system( MovCmd ) 
