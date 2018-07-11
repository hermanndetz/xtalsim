/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/

#ifndef __SIMULATION_BOX_RENDERER__
#define __SIMULATION_BOX_RENDERER__

#include "projectConfigure.h"

#ifdef __VTK__

#include <memory>
#include <sstream>
#include <algorithm>
#include <vector>
#include <tuple>
#include <iostream>
#include <unordered_map>
#include <cctype>
#include <string>

#include <vtkVersion.h>
#include <vtkConfigure.h>
#include <vtkObjectFactory.h>
#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolygon.h>
#include <vtkPolyLine.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkProperty.h>
#include <vtkFloatArray.h>
#include <vtkLookupTable.h>
#include <vtkGlyph3D.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSphereSource.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkAxesActor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkLine.h>
#include <vtkCamera.h>
#include <vtkLegendBoxActor.h>
#include <vtkNamedColors.h>
#include <vtkColorTransferFunction.h>
#include <vtkPlaneSource.h>
#include <vtkWindowToImageFilter.h>
#include <vtkJPEGWriter.h>
#include <vtkPNGWriter.h>
#include <vtkCommand.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkActorCollection.h>
#include <vtkRendererCollection.h>

#include <VtkDefs.h>

// For compatibility with new VTK generic data arrays
#ifdef vtkGenericDataArray_h
#define InsertNextTupleValue InsertNextTypedTuple
#endif

#include <misc/Configuration.h>
#include <physics/Field3D.h>
#include <physics/MaterialCollection.h>
#include <physics/PeriodicTable.h>
#include <physics/SimulationBox.h>
#include <physics/Atom.h>
#include <physics/Vector3D.h>

//##############################################################################

class SimulationBoxRenderer {
private:
    const std::string logName_;

    const uint16_t windowWidth_;      //!< Width of VTK window in pixel.
    const uint16_t windowHeight_;     //!< Height of VTK window in pixel.
    
    vtkSmartPointer<vtkRenderer> renderer_;
    vtkSmartPointer<vtkRenderWindow> renderWindow_;
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor_;

    Configuration config_;

    SimulationBox *simbox_;
    const MaterialCollection *materials_;
    
    //! stores intial camera positions for each renderer
    std::unordered_map<vtkRenderer*, Vector3D<spaceType>> initialCameraPosition_;
    std::unordered_map<vtkRenderer*, Vector3D<spaceType>> initialCameraViewUp_;
    std::unordered_map<vtkRenderer*, Range3D<indexType>> displayRange_;

    //! Bounds of renderer, which was called last
    double bounds_[6]{0.0};

    //! Information if particular layer contains modified atoms
    std::vector<bool> modifiedAtomsInLayer_{};
    //! Information if a layer contains only modified atoms
    std::vector<bool> modifiedAtomsOnlyInLayer_{};

    //! Determines, which layers (exclusively) contain modified atoms
    void evaluateModifiedAtomsInLayers (uint32_t modificationIndexLimit=0);
    //! Checks, if mofified flags were evaluated already
    bool modifiedAtomsSet (void) const;

    void addCrystalViews(double commandLineHeightShare,
                         const std::shared_ptr<Field3D<double>> strainField);
    void addInterfaces(double commandLineHeightShare,
                       const std::vector<Range3D<indexType>> interfaces);

    //! Determines if atom shall be processed under current settings.
    inline bool continueLoop(Atom *atom) const;
    
public:
    //! Constructor
    SimulationBoxRenderer (const Configuration &config, const char *logName="SimulationBoxRenderer");

    //! Destructor
    ~SimulationBoxRenderer ();

    void clear (void);
    void render (SimulationBox &simbox);

    void renderAxes3D ();
    void renderBoundingBox (const SimulationBox &simbox);
    void start ();
    void renderQuadViewports (SimulationBox &simbox,
                              const MaterialCollection &materials,
                              const Configuration &config);


    void renderAtoms (vtkSmartPointer<vtkRenderer> renderer);
    void renderAtomsWithinLimit(vtkRenderer *renderer, const int limit);
    
    void renderBonds (vtkSmartPointer<vtkRenderer> renderer,
                      const bool strain=false,
                      const double refLatticeConst=-1.0);
    void renderAbsStrainPlanes (vtkSmartPointer<vtkRenderer> renderer,
                                const SimulationBox &simbox,
                                const std::shared_ptr<Field3D<double>> strainField) const;

    void renderPerspectiveWithInterfaces (SimulationBox &simbox,
                                          const std::vector<Range3D<indexType>> interfaces,
                                          const MaterialCollection &materials,
                                          const Configuration &config,
                                          std::shared_ptr<Field3D<double>> strainField = nullptr,
                                          bool persist=true);


    void resetCamera(vtkRenderer *renderer);
    uint32_t getMaxModificationIndex(void);
    void saveImage (const std::string filename) const;

};

//##############################################################################

//class CommandInteractorStyle;

class KeyPressInteractorStyle : public vtkInteractorStyleTrackballCamera
{
private:
    std::string logName_ = "KeyPressInteractor";
    
    //! Width and height of interaction Window.
    uint16_t windowWidth_, windowHeight_;

    SimulationBoxRenderer *simboxRenderer_;
    vtkTextActor *commandLineActor_;

    uint32_t animationTimeout_=200;
    std::string imagePreamble_="";
    
    std::string commandString_;
    void updateCommandLine(void);

    bool commandLineMode_ = false;

    //! toggle display of bonds in all renderers
    void toggleBonds(vtkRenderWindowInteractor *rwi);
    void animate(vtkRenderWindowInteractor *rwi);
    
public:
    static KeyPressInteractorStyle* New();
    vtkTypeMacro(KeyPressInteractorStyle, vtkInteractorStyleTrackballCamera);

    //void SetCommandInteractor(CommandInteractorStyle *commandInteractor);
    
    void setSize(const uint16_t windowWidth, const uint16_t windowHeight);
    void setSimboxRenderer(SimulationBoxRenderer *simboxRenderer);
    void setCommandLine(vtkTextActor *textActor);
    void setLogName(std::string logName);
    void runCommand(const std::string command);
    
    void OnKeyPress() VTK_OVERRIDE;
    void OnChar() VTK_OVERRIDE;
    void OnLeftButtonDown() VTK_OVERRIDE;
    /*  void OnMouseMove() VTK_OVERRIDE;
      void OnLeftButtonDown() VTK_OVERRIDE;
      void OnLeftButtonUp() VTK_OVERRIDE;
      void OnMiddleButtonDown() VTK_OVERRIDE;
      void OnMiddleButtonUp() VTK_OVERRIDE;
      void OnRightButtonDown() VTK_OVERRIDE;
      void OnRightButtonUp() VTK_OVERRIDE;
      void OnMouseWheelForward() VTK_OVERRIDE;
      void OnMouseWheelBackward() VTK_OVERRIDE;*/
};

//##############################################################################

#endif

#endif

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
