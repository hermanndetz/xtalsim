/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/

#include "SimulationBoxRenderer.h"

#ifdef __VTK__

SimulationBoxRenderer::SimulationBoxRenderer (const Configuration &config, const char *logName):
    config_(config), windowWidth_(1280), windowHeight_(864)
{

    renderer_ = vtkSmartPointer<vtkRenderer>::New();
    renderWindow_ = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindowInteractor_ = vtkSmartPointer<vtkRenderWindowInteractor>::New();

    // create a renderer, render window, and interactor
    renderWindow_->AddRenderer(renderer_);
    renderWindowInteractor_->SetRenderWindow(renderWindow_);

    renderer_->SetBackground(255, 255, 255);
    renderWindow_->SetSize(864,864);

    // render and interact
    vtkSmartPointer<KeyPressInteractorStyle> style =
        vtkSmartPointer<KeyPressInteractorStyle>::New(); //like paraview

    style->setSize(windowWidth_, windowHeight_);
    style->setSimboxRenderer(this);

    renderWindowInteractor_->SetInteractorStyle( style );
}

//------------------------------------------------------------------------------

SimulationBoxRenderer::~SimulationBoxRenderer () {

}

//------------------------------------------------------------------------------

void SimulationBoxRenderer::clear(void) {
    renderer_->RemoveAllViewProps();
}

//------------------------------------------------------------------------------

void SimulationBoxRenderer::render (SimulationBox &simbox) 
{
    simbox_ = &simbox;
    
    renderAtoms(renderer_);
}

//******************************************************************************

//! \param atom Atom that shall be checked.
//! \return true if atom shall be skipped
inline bool SimulationBoxRenderer::continueLoop(Atom *atom) const
{

    uint8_t outOfPlaneDimension = simbox_->getOutOfPlaneDimension();
    
    // contains all abort conditions for a single atom
    if (modifiedAtomsInLayer_[atom->getIndex()[outOfPlaneDimension]] == false) {
        return true;
    }
        
    // required to deal with old XML files
    if ((config_.modifiedOnly == true) && (config_.filterModificationState.size() == 0)) {
        if (modifiedAtomsOnlyInLayer_[atom->getIndex()[outOfPlaneDimension]])
            return true;

        if (atom->wasModified() == config_.modifiedNegative)
            return true;

        if (atom->wasModified() &&
            (config_.modificationIndexMax > 0) &&
            (atom->getModificationOrder() > config_.modificationIndexMax)){

            return true;
        }
    }

    if (config_.filterModificationState.size() > 0) {
        std::bitset<AtomState::StateCount> tmpFlags = atom->getState() &
            config_.getModificationState();


        if (tmpFlags.none())
            return true;
    }
    
    return false;
}
    
//------------------------------------------------------------------------------

// void SimulationBoxRenderer::addAtomsModificationIndex (vtkSmartPointer<vtkRenderer> renderer,
//                                                        int minIndex, int maxIndex)
// {

//     std::vector<Element>elements = periodicTable_->get();
//     vtkSmartPointer<vtkLookupTable> colors = vtkSmartPointer<vtkLookupTable>::New();
//     std::vector<std::tuple<uint8_t,elementType>> colorElementRelation;

//     for (int i = 0; i < elements.size(); i++) {
//         colors->SetTableValue(i, elements[i].color.get(0),
//                               elements[i].color.get(1), elements[i].color.get(2), 1.0);
//         colorElementRelation.push_back(std::make_tuple(i, elements[i].id));
//     }
    
//     for (auto atom: simbox_->getLattice().getAtomList(displayRange_[renderer])) {
//         if ( continueLoop(atom))
//             continue;

//         Vector3D<spaceType> pos = atom->getPosition();
//         points->InsertNextPoint(pos[0], pos[1], pos[2]);
//         scales->InsertNextValue(0.7); // fixed radius for now

//         auto colElement = std::find_if(colorElementRelation.begin(), colorElementRelation.end(), 
//                                        [=](const std::tuple<uint8_t,elementType>& val) {
//                                            return std::get<1>(val) == atom->getElementId();});

//         auto colID = std::distance(colorElementRelation.begin(), colElement);

//         col->InsertNextValue(colID);
//     }
    
// }

//------------------------------------------------------------------------------

void SimulationBoxRenderer::renderAtoms (vtkSmartPointer<vtkRenderer> renderer)
{
    uint8_t outOfPlaneDimension = simbox_->getOutOfPlaneDimension();
    std::vector<Element>elements = PeriodicTable::getInstance().get();
    std::vector<std::tuple<uint8_t,elementType>> colorElementRelation;
   
    // create points
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    // setup scales
    vtkSmartPointer<vtkFloatArray> scales = vtkSmartPointer<vtkFloatArray>::New();
    scales->SetName("scales");

     // setup color label
    vtkSmartPointer<vtkFloatArray> col = vtkSmartPointer<vtkFloatArray>::New();
    col->SetName("col");

    vtkSmartPointer<vtkLookupTable> colors = vtkSmartPointer<vtkLookupTable>::New();
    // 4 table values: (index), red, green, blue, opacity
    colors->SetNumberOfTableValues(elements.size());

    for (int i = 0; i < elements.size(); i++) {
        colors->SetTableValue(i, elements[i].color.get(0),
                              elements[i].color.get(1), elements[i].color.get(2), 1.0);
        colorElementRelation.push_back(std::make_tuple(i, elements[i].id));
    }

    if (modifiedAtomsSet() == false)
        evaluateModifiedAtomsInLayers();

    for (auto atom: simbox_->getLattice().getAtomList(displayRange_[renderer])) {
        if ( continueLoop(atom) )
            continue;

        Vector3D<spaceType> pos = atom->getPosition();
        points->InsertNextPoint(pos[0], pos[1], pos[2]);
        scales->InsertNextValue(0.7); // fixed radius for now

        auto colElement = std::find_if(colorElementRelation.begin(), colorElementRelation.end(), 
                                       [=](const std::tuple<uint8_t,elementType>& val) {
                                           return std::get<1>(val) == atom->getElementId();});

        auto colID = std::distance(colorElementRelation.begin(), colElement);

        col->InsertNextValue(colID);
    }

    // grid structured to append center, radius and color label
    vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    grid->SetPoints(points);
    grid->GetPointData()->AddArray(scales);
    grid->GetPointData()->SetActiveScalars("scales"); // !!!to set radius first
    grid->GetPointData()->AddArray(col);

    grid->GetBounds(bounds_);

    // create anything you want here, we will use a sphere for the demo
    vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
    sphereSource->SetThetaResolution(16);
    sphereSource->SetPhiResolution(16);

    // object to group sphere and grid and keep smooth interaction
    vtkSmartPointer<vtkGlyph3D> glyph3D = vtkSmartPointer<vtkGlyph3D>::New();
    glyph3D->SetInputData(grid);
    glyph3D->SetSourceConnection(sphereSource->GetOutputPort());

    // create a mapper and actor
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(glyph3D->GetOutputPort());

    mapper->SetScalarModeToUsePointFieldData(); // without, color are displayed regarding radius and not color label
    mapper->SetScalarRange(0, elements.size()); // to scale color label (without, col should be between 0 and 1)
    mapper->SelectColorArray("col"); // !!!to set color (nevertheless you will have nothing)
    mapper->SetLookupTable(colors);

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    renderer->AddActor(actor);

}

//------------------------------------------------------------------------------

void SimulationBoxRenderer::renderBonds (vtkSmartPointer<vtkRenderer> renderer,
                                         const bool strain,
                                         const double refLatticeConst) {

    const Range3D<indexType> &range = displayRange_[renderer];
    std::vector<std::tuple<Vector3D<indexType>, Vector3D<indexType>, double>> strainInfo;
    uint8_t outOfPlaneDimension = simbox_->getOutOfPlaneDimension();
    Lattice &lattice = simbox_->getLattice();
   
    if (strain == true) {
        simbox_->generateNeighbors(1, 0.0, false);
        strainInfo = simbox_->getStrainInfo(*materials_, refLatticeConst);
    }

    // prepare color gradient
    vtkSmartPointer<vtkColorTransferFunction> ctf = vtkSmartPointer<vtkColorTransferFunction>::New();
    ctf->AddRGBPoint(-0.1, (double)215/(double)255, (double)48/(double)255, (double)39/(double)255);
    ctf->AddRGBPoint(0.0, 0.0, 0.0, 0.0);
    ctf->AddRGBPoint(0.1, (double)69/(double)255, (double)117/(double)255, (double)180/(double)255);

    // create points and lines
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();

    colors->SetNumberOfComponents(3);

    if (modifiedAtomsSet() == false)
        evaluateModifiedAtomsInLayers();
   
    uint32_t i=0;
    for (auto atom: lattice.getAtomList(range)) {

        if ( continueLoop(atom) )
            continue;
        
        Vector3D<spaceType> startPos = atom->getPosition();

        for (auto neighbor: atom->getNeighbors()) {
            Atom *neighborAtom = lattice(neighbor);

            if ( ! neighbor.isInRange(range) )
                continue;

            if ( continueLoop(neighborAtom) )
                continue;

            Vector3D<spaceType> endPos = neighborAtom->getPosition();

            points->InsertNextPoint(startPos[0], startPos[1], startPos[2]);
            points->InsertNextPoint(endPos[0], endPos[1], endPos[2]);

            unsigned char bondCol[3] = {0,0,0};

            if (strain == true) {
                double tmpStrain = 0.0;
                auto tmpVectorTuple = std::make_tuple(atom->getIndex(), neighbor, 0.0);

                std::vector<std::tuple<Vector3D<indexType>, Vector3D<indexType>, double>>::const_iterator  bondStrain =
                    std::find_if(strainInfo.begin(), strainInfo.end(),
                                 [&tmpVectorTuple](const std::tuple<Vector3D<indexType>, Vector3D<indexType>, double>& item) -> 
                                 bool {
                                     if (((std::get<0>(item) == std::get<0>(tmpVectorTuple)) &&
                                          (std::get<1>(item) == std::get<1>(tmpVectorTuple))) ||
                                         ((std::get<0>(item) == std::get<1>(tmpVectorTuple)) &&
                                          (std::get<1>(item) == std::get<0>(tmpVectorTuple)))) {
                                         return true;
                                     }
                                     return false;});

                if (bondStrain != strainInfo.end()) {
                    tmpStrain = std::get<2>(*bondStrain);
                }

                double bondColor[3] = {0.0, 0.0, 0.0};
                ctf->GetColor(tmpStrain, bondColor);

                bondCol[0] = (unsigned char)(bondColor[0]*255.0);
                bondCol[1] = (unsigned char)(bondColor[1]*255.0);
                bondCol[2] = (unsigned char)(bondColor[2]*255.0);

            }

            colors->InsertNextTupleValue(bondCol);
            
            vtkSmartPointer<vtkLine> bond = vtkSmartPointer<vtkLine>::New();
            bond->GetPointIds()->SetId(0, i++);
            bond->GetPointIds()->SetId(1, i++);
            lines->InsertNextCell(bond);
        }

    }
    
    vtkSmartPointer<vtkPolyData> bonds = vtkSmartPointer<vtkPolyData>::New();
    bonds->SetPoints(points);
    bonds->SetLines(lines);
    bonds->GetCellData()->SetScalars(colors);

    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(bonds);

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetLineWidth(2);
    renderer->AddActor(actor);
    
}

//------------------------------------------------------------------------------

void SimulationBoxRenderer::renderAbsStrainPlanes (vtkSmartPointer<vtkRenderer> renderer,
                                                   const SimulationBox &simbox,
                                                   const std::shared_ptr<Field3D<double>> strainField) const {

    Vector3D<spaceType> latticeSize = simbox.getSize();
    uint32_t downsamplingFactor[3] = {4, 4, 4};

    Field3D<double> strainFieldDS = strainField->downsampling(downsamplingFactor[0],
                                                              downsamplingFactor[1], downsamplingFactor[2],
                                                              DownsamplingFieldOperation::Avg);

    // render plane by plane
    // (stop either when there are no more strain data or atomic planes)
    for (uint32_t k = 0; k < simbox.getLattice().getSize()[2]; k+=4) {

        vtkSmartPointer<vtkPlaneSource> currentPlane = vtkSmartPointer<vtkPlaneSource>::New();

        double avgPlanePosition = simbox.getLattice().getAvgCoordinate(k, 2, config_.modifiedOnly);

        // in perfect structures, this will also filter out plane 0
        // \todo needs a better evaluation in future.
        if (std::fpclassify(avgPlanePosition) == FP_ZERO)
            continue;

        currentPlane->SetOrigin(0.0, 0.0, avgPlanePosition);
        currentPlane->SetPoint1(latticeSize[0], 0.0, avgPlanePosition);
        currentPlane->SetPoint2(0.0, latticeSize[1], avgPlanePosition);
        currentPlane->SetXResolution((simbox.getLattice().getSize()[0]-1)/downsamplingFactor[0]);
        currentPlane->SetYResolution((simbox.getLattice().getSize()[1]-1)/downsamplingFactor[1]);
        currentPlane->Update();

        // prepare color gradient
        double colorCenter = (strainField->max() - strainField->min())/2.0;

        if ((strainField->max() > 0.0) && (strainField->min() < 0.0))
            colorCenter = 0.0;


        vtkSmartPointer<vtkUnsignedCharArray> cellData = vtkSmartPointer<vtkUnsignedCharArray>::New();
        cellData->SetNumberOfComponents(3);
        cellData->SetNumberOfTuples(currentPlane->GetOutput()->GetNumberOfCells());

        Vector3D<indexType> fieldSize = strainField->getSize();

        vtkSmartPointer<vtkColorTransferFunction> ctf = vtkSmartPointer<vtkColorTransferFunction>::New();
        ctf->AddRGBPoint(strainFieldDS.min(), (double)215/(double)255, (double)48/(double)255, (double)39/(double)255);
        //ctf->AddRGBPoint(-10, (double)215/(double)255, (double)48/(double)255, (double)39/(double)255);
        ctf->AddRGBPoint(colorCenter, 0.0, 0.0, 0.0);
        ctf->AddRGBPoint(strainFieldDS.max(), (double)69/(double)255, (double)117/(double)255, (double)180/(double)255);
        //ctf->AddRGBPoint(10, (double)69/(double)255, (double)117/(double)255, (double)180/(double)255);

        double rgb[3];

        uint32_t l = 0;

        Vector3D<indexType> dsSize = strainFieldDS.getSize();

        for (uint32_t j = 0; j < dsSize[1]; j+=1) {
            for (uint32_t i = 0; i < dsSize[0]; i+=1) {
                ctf->GetColor(strainFieldDS(i, j, k/downsamplingFactor[2]), rgb); 

                rgb[0] = rgb[0] * 255.0;
                rgb[1] = rgb[1] * 255.0;
                rgb[2] = rgb[2] * 255.0;

                cellData->InsertTuple(l++, rgb);
            }
        }

        currentPlane->GetOutput()->GetCellData()->SetScalars(cellData);

        vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputConnection(currentPlane->GetOutputPort());

        vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);
        actor->GetProperty()->SetOpacity(0.67);

        renderer->AddActor(actor);
    }
}

//******************************************************************************

void SimulationBoxRenderer::renderAxes3D () {
  // Setup four points
  vtkSmartPointer<vtkPoints> points =
    vtkSmartPointer<vtkPoints>::New();
  points->InsertNextPoint(-2.0, -2.0, -2.0);
  points->InsertNextPoint(1.0, -2.0, -2.0);
  points->InsertNextPoint(-2.0, 1.0, -2.0);
  points->InsertNextPoint(-2.0, -2.0, 1.0);

  // Create the polygon
  vtkSmartPointer<vtkPolyLine> polygon =
    vtkSmartPointer<vtkPolyLine>::New();
  polygon->GetPointIds()->SetNumberOfIds(7); //make a quad
  polygon->GetPointIds()->SetId(0, 0);
  polygon->GetPointIds()->SetId(1, 1);
  polygon->GetPointIds()->SetId(2, 0);
  polygon->GetPointIds()->SetId(3, 2);
  polygon->GetPointIds()->SetId(4, 0);
  polygon->GetPointIds()->SetId(5, 3);
  polygon->GetPointIds()->SetId(6, 0);

  // Add the polygon to a list of polygons
  vtkSmartPointer<vtkCellArray> polygons =
    vtkSmartPointer<vtkCellArray>::New();
  polygons->InsertNextCell(polygon);

  unsigned char colorDef[3] = {0, 0, 0};

  vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
  colors->SetNumberOfComponents(3);
  colors->SetName("Colors");
  colors->InsertNextTupleValue(colorDef);
  colors->InsertNextTupleValue(colorDef);
  colors->InsertNextTupleValue(colorDef);
  colors->InsertNextTupleValue(colorDef);
  colors->InsertNextTupleValue(colorDef);
  colors->InsertNextTupleValue(colorDef);
  colors->InsertNextTupleValue(colorDef);

  // Create a PolyData
  vtkSmartPointer<vtkPolyData> polygonPolyData =
    vtkSmartPointer<vtkPolyData>::New();
  polygonPolyData->SetPoints(points);
  polygonPolyData->SetLines(polygons);
  //polygonPolyData->SetPolys(polygons);
  polygonPolyData->GetPointData()->SetScalars(colors);

  // Create a mapper and actor
  vtkSmartPointer<vtkPolyDataMapper> mapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
  mapper->SetInput(polygonPolyData);
#else
  mapper->SetInputData(polygonPolyData);
#endif

  vtkSmartPointer<vtkActor> actor =
    vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);

  // Visualize

  renderer_->AddActor(actor);
}

//------------------------------------------------------------------------------

void SimulationBoxRenderer::renderBoundingBox (const SimulationBox &simbox) {
    Vector3D<double> origin = {0.0, 0.0, 0.0};
    Vector3D<double> spaceDiagonal = simbox.getSize();

    int numberOfPolygonPoints{0};

  // Setup four points
  vtkSmartPointer<vtkPoints> points =
    vtkSmartPointer<vtkPoints>::New();
  points->InsertNextPoint(origin[0], origin[1], origin[2]);
  points->InsertNextPoint(spaceDiagonal[0], origin[1], origin[2]);
  points->InsertNextPoint(spaceDiagonal[0], spaceDiagonal[1], origin[2]);
  points->InsertNextPoint(origin[0], spaceDiagonal[1], origin[0]);
  points->InsertNextPoint(origin[0], origin[1], spaceDiagonal[2]);
  points->InsertNextPoint(spaceDiagonal[0], origin[1], spaceDiagonal[2]);
  points->InsertNextPoint(spaceDiagonal[0], spaceDiagonal[1], spaceDiagonal[2]);
  points->InsertNextPoint(origin[0], spaceDiagonal[1], spaceDiagonal[2]);

  vtkSmartPointer<vtkCellArray> polygons =
    vtkSmartPointer<vtkCellArray>::New();

  // Create the base polygon
  vtkSmartPointer<vtkPolyLine> polyBase =
    vtkSmartPointer<vtkPolyLine>::New();
  polyBase->GetPointIds()->SetNumberOfIds(5); //make a quad
  polyBase->GetPointIds()->SetId(0, 0);
  polyBase->GetPointIds()->SetId(1, 1);
  polyBase->GetPointIds()->SetId(2, 2);
  polyBase->GetPointIds()->SetId(3, 3);
  polyBase->GetPointIds()->SetId(4, 0);
  numberOfPolygonPoints += 5;

  polygons->InsertNextCell(polyBase);


  // Create the top polygon
  vtkSmartPointer<vtkPolyLine> polyTop =
    vtkSmartPointer<vtkPolyLine>::New();
  polyTop->GetPointIds()->SetNumberOfIds(5); //make a quad
  polyTop->GetPointIds()->SetId(0, 4);
  polyTop->GetPointIds()->SetId(1, 5);
  polyTop->GetPointIds()->SetId(2, 6);
  polyTop->GetPointIds()->SetId(3, 7);
  polyTop->GetPointIds()->SetId(4, 4);
  numberOfPolygonPoints += 5;

  polygons->InsertNextCell(polyTop);

  // Create the sidewalls polygon
  vtkSmartPointer<vtkLine> line[4] = {vtkSmartPointer<vtkLine>::New()};


  for (int i: {0, 1, 2, 3})
  {
      line[i]->GetPointIds()->SetId(0,i);
      line[i]->GetPointIds()->SetId(1,i+4);
      numberOfPolygonPoints += 2;
      polygons->InsertNextCell(line[i]);
  }


  unsigned char colorDef[3] = {0, 0, 0};

  vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
  colors->SetNumberOfComponents(3);
  colors->SetName("Colors");

  for (int i = 0; i < numberOfPolygonPoints; i++) {
      colors->InsertNextTupleValue(colorDef);
  }

  // Create a PolyData
  vtkSmartPointer<vtkPolyData> polygonPolyData =
    vtkSmartPointer<vtkPolyData>::New();
  polygonPolyData->SetPoints(points);
  polygonPolyData->SetLines(polygons);
  //polygonPolyData->SetPolys(polygons);
  polygonPolyData->GetPointData()->SetScalars(colors);

  // Create a mapper and actor
  vtkSmartPointer<vtkPolyDataMapper> mapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
  mapper->SetInput(polygonPolyData);
#else
  mapper->SetInputData(polygonPolyData);
#endif

  vtkSmartPointer<vtkActor> actor =
    vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);

  // Visualize

  renderer_->AddActor(actor);
}

//------------------------------------------------------------------------------

void SimulationBoxRenderer::start ()
{
    renderWindow_->Render();
    renderWindowInteractor_->Start();
}

//------------------------------------------------------------------------------

void SimulationBoxRenderer::renderQuadViewports (SimulationBox &simbox,
                                                 const MaterialCollection &materials,
                                                 const Configuration &config) {
    config_ = config;
    simbox_ = &simbox;
    materials_ = &materials;

    auto largestExtent = [](double *boundsArray) {
        double result = 0.0;

        for (auto i: {0, 1, 2}) {
            if ((boundsArray[i*2+1] - boundsArray[i*2]) > result)
                result = (boundsArray[i*2+1] - boundsArray[i*2]);
        }

        return result;
    };

    double xmins[4] = {0.0, 0.0, 0.0, 0.25};
    double xmaxs[4] = {0.25, 0.25, 0.25, 1.0};
    double ymins[4] = {0.0, 0.33, 0.66, 0.0};
    double ymaxs[4] = {0.33, 0.66, 1.0, 1.0};

    Vector3D<double> origin = {0.0, 0.0, 0.0};
    Vector3D<double> spaceDiagonal = simbox.getSize();
    Vector3D<double> cameraPositions[4] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0},
                                           {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

    cameraPositions[0][0] = (spaceDiagonal[0] - origin[0]) / 2.0;
    cameraPositions[0][1] = (spaceDiagonal[1] - origin[1]) / 2.0;
    cameraPositions[0][2] = (spaceDiagonal[2] - origin[2]) * -5.0;

    cameraPositions[1][0] = (spaceDiagonal[0] - origin[0]) / 2.0;
    cameraPositions[1][1] = (spaceDiagonal[1] - origin[1]) * -5.0;
    cameraPositions[1][2] = (spaceDiagonal[2] - origin[2]) / 2.0;

    cameraPositions[2][0] = (spaceDiagonal[0] - origin[0]) * -5.0;
    cameraPositions[2][1] = (spaceDiagonal[1] - origin[1]) / 2.0;
    cameraPositions[2][2] = (spaceDiagonal[2] - origin[2]) / 2.0;

    cameraPositions[3][0] = (spaceDiagonal[0] - origin[0]) * -5.5;
    cameraPositions[3][1] = (spaceDiagonal[1] - origin[1]) * -4.5;
    cameraPositions[3][2] = (spaceDiagonal[2] - origin[2]) * 1.30;

    for (int i: {0, 1, 2, 3}) {
        vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();

        renderWindow_->AddRenderer(renderer);
        renderer->SetViewport(xmins[i], ymins[i], xmaxs[i], ymaxs[i]);
        renderer->SetBackground(1.0, 1.0, 1.0);

        vtkSmartPointer<vtkCamera> camera = vtkSmartPointer<vtkCamera>::New();
        camera->SetPosition(cameraPositions[i][0], cameraPositions[i][1], cameraPositions[i][2]);

        vtkSmartPointer<vtkLegendBoxActor> legend = vtkSmartPointer<vtkLegendBoxActor>::New();
        vtkSmartPointer<vtkNamedColors> colors = vtkSmartPointer<vtkNamedColors>::New();
        double color[4];

        if (i < 3) {
            camera->ParallelProjectionOn();
            renderer->InteractiveOff();

            legend->SetNumberOfEntries(1);

            legend->GetPositionCoordinate()->SetCoordinateSystemToView();
            legend->GetPositionCoordinate()->SetValue(0.5, -1.0);

            legend->GetPosition2Coordinate()->SetCoordinateSystemToView();
            legend->GetPosition2Coordinate()->SetValue(1.0, -0.5);

            legend->UseBackgroundOff();
            legend->BorderOff();

            switch (i) {
            case 0:
                camera->SetFocalPoint(cameraPositions[i][0], cameraPositions[i][1], 0.0);
                camera->SetParallelScale(10.0);
                renderer->SetBackground(1.0, 1.0, 1.0);
                camera->SetViewUp(-1.0, 0.0, 0.0);

                colors->GetColor("black", color);
                legend->SetEntryString(0, "(0 0 1)");
                legend->SetEntryColor(0, color[0], color[1], color[2]);
                break;
            case 1:
                renderer->SetBackground(0.8, 0.8, 0.8);

                camera->SetFocalPoint(cameraPositions[i][0], 0.0, cameraPositions[i][2]);
                camera->SetParallelScale(10.0);
                camera->ParallelProjectionOn();
                camera->SetViewUp(0.0, 0.0, 1.0);
                camera->Azimuth(1.0);
                camera->SetViewUp(0.0,-1.0,0.0);

                colors->GetColor("black", color);
                legend->SetEntryString(0, "(0 1 0)");
                legend->SetEntryColor(0, color[0], color[1], color[2]);

                break;
            case 2:
                camera->SetFocalPoint(0.0, cameraPositions[i][1], cameraPositions[i][2]);
                camera->SetParallelScale(10.0);
                camera->ParallelProjectionOn();
                renderer->SetBackground(1.0, 1.0, 1.0);
                camera->SetViewUp(0.0,1.0,0.0);

                colors->GetColor("black", color);
                legend->SetEntryString(0, "(1 0 0)");
                legend->SetEntryColor(0, color[0], color[1], color[2]);
                break;
            }


            renderer->AddActor(legend);
        } else {
            camera->SetFocalPoint(spaceDiagonal[0]/2.0, spaceDiagonal[1]/2.0, spaceDiagonal[2]/2.0);
            camera->SetViewUp(0.0, 0.0, 1.0);
        }

        renderer->SetActiveCamera(camera);
        renderer->ResetCamera();
        displayRange_[renderer] = Range3D<indexType>();

        if (config.renderBonds == true) {
            if (i == 3)
                renderBonds(renderer, config.visualizeBondStrain, config.refLatticeConstant);
            else
                renderBonds(renderer);
        }

        renderAtoms(renderer);

        auto le = largestExtent(bounds_);

        // autoscale individual views to fill the
        // individual render windows
        // camera 3 gets auto scale through ResetCamera function.
        if (i<3)  camera->SetParallelScale(1.0*le);

    }

    renderWindow_->Render();
    renderWindowInteractor_->Start();

    config_ = Configuration();
}

//******************************************************************************

void SimulationBoxRenderer::addCrystalViews(double commandLineHeightShare,
                                            const std::shared_ptr<Field3D<double>> strainField)
{
    double xmins[4] = {0.0, 0.0, 0.0, 0.25};
    double xmaxs[4] = {0.25, 0.25, 0.25, 0.75};
    double ymins[4] = {commandLineHeightShare,
                       commandLineHeightShare + (1.0-commandLineHeightShare)/3,
                       commandLineHeightShare + (1.0-commandLineHeightShare)*2/3,
                       commandLineHeightShare};
    double ymaxs[4] = {commandLineHeightShare + (1.0-commandLineHeightShare)/3,
                       commandLineHeightShare + (1.0-commandLineHeightShare)*2/3, 1.0, 1.0};

    Vector3D<spaceType> origin = {0.0, 0.0, 0.0};
    Vector3D<spaceType> spaceDiagonal = simbox_->getSize();
    Vector3D<spaceType> cameraPositions[4] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0},
                                              {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

    auto largestExtent = [](double *boundsArray) {
        double result = 0.0;

        for (auto i: {0, 1, 2}) {
            if ((boundsArray[i*2+1] - boundsArray[i*2]) > result)
                result = (boundsArray[i*2+1] - boundsArray[i*2]);
        }

        return result;
    };

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    for (int i: {0, 1, 2}){

        for (int j: {0, 1, 2})
            cameraPositions[i][j] = (spaceDiagonal[j] - origin[j]) / 2.0;

        cameraPositions[i][2-i] *= -10.0;
    }
    
    for (int j: {0, 1, 2}){
        cameraPositions[3][j] = (spaceDiagonal[j] - origin[j]) / 2;
        cameraPositions[3][j] += spaceDiagonal[0] * 4 *
            config_.cameraDirection[j] * (1/config_.cameraZoom);
    }
        
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    for (int i: {0, 1, 2, 3}) {
        Vector3D<spaceType> tmpPosition;
        vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();

        renderWindow_->AddRenderer(renderer);
        renderer->SetViewport(xmins[i], ymins[i], xmaxs[i], ymaxs[i]);
        renderer->SetBackground(1.0, 1.0, 1.0);

        vtkSmartPointer<vtkCamera> camera = vtkSmartPointer<vtkCamera>::New();
        tmpPosition  = {cameraPositions[i][0], cameraPositions[i][1], cameraPositions[i][2]};
        camera->SetPosition(tmpPosition[0], tmpPosition[1], tmpPosition[2]);
        initialCameraPosition_[renderer] = tmpPosition;

        vtkSmartPointer<vtkLegendBoxActor> legend = vtkSmartPointer<vtkLegendBoxActor>::New();
        vtkSmartPointer<vtkNamedColors> colors = vtkSmartPointer<vtkNamedColors>::New();
        double color[4];

        if (i < 3) {
            camera->ParallelProjectionOn();
            renderer->InteractiveOff();

            legend->SetNumberOfEntries(1);

            legend->GetPositionCoordinate()->SetCoordinateSystemToView();
            legend->GetPositionCoordinate()->SetValue(0.5, -1.0);

            legend->GetPosition2Coordinate()->SetCoordinateSystemToView();
            legend->GetPosition2Coordinate()->SetValue(1.0, -0.5);

            legend->UseBackgroundOff();
            legend->BorderOff();

            switch (i) {
            case 0:
                camera->SetFocalPoint(cameraPositions[i][0], cameraPositions[i][1], 0.0);
                camera->SetParallelScale(10.0);
                renderer->SetBackground(1.0, 1.0, 1.0);

                tmpPosition = {-1.0, 0.0, 0.0};
                camera->SetViewUp(tmpPosition[0], tmpPosition[1], tmpPosition[2]);
                initialCameraViewUp_[renderer] = tmpPosition;
                
                colors->GetColor("black", color);
                legend->SetEntry(0, (vtkPolyData *)NULL, "(0 0 1)", color);

                break;
            case 1:
                renderer->SetBackground(0.8, 0.8, 0.8);

                camera->SetFocalPoint(cameraPositions[i][0], 0.0, cameraPositions[i][2]);
                camera->SetParallelScale(10.0);
                camera->ParallelProjectionOn();
                camera->SetViewUp(0.0, 0.0, 1.0);
                camera->Azimuth(1.0);

                tmpPosition = {0.0, -1.0, 0.0};
                camera->SetViewUp(tmpPosition[0], tmpPosition[1], tmpPosition[2]);
                initialCameraViewUp_[renderer] = tmpPosition;
                
                colors->GetColor("black", color);
                legend->SetEntry(0, (vtkPolyData *)NULL, "(0 1 0)", color);

                break;
            case 2:
                camera->SetFocalPoint(0.0, cameraPositions[i][1], cameraPositions[i][2]);
                camera->SetParallelScale(10.0);
                camera->ParallelProjectionOn();
                renderer->SetBackground(1.0, 1.0, 1.0);

                tmpPosition = {0.0, 1.0, 0.0};
                camera->SetViewUp(tmpPosition[0], tmpPosition[1], tmpPosition[2]);
                initialCameraViewUp_[renderer] = tmpPosition;
                
                colors->GetColor("black", color);
                legend->SetEntry(0, (vtkPolyData *)NULL, "(1 0 0)", color);

                break;
            }


            renderer->AddActor(legend);
        } else {
            if ( config_.cameraFocalPoint == 0 )
                camera->SetFocalPoint(spaceDiagonal[0]/2.0, spaceDiagonal[1]/2.0, spaceDiagonal[2]/2.0);
            else
                camera->SetFocalPoint(config_.cameraFocalPoint[0], config_.cameraFocalPoint[1],
                                      config_.cameraFocalPoint[2]);

            tmpPosition = {0.0, 0.0, 1.0};
            camera->SetViewUp(tmpPosition[0], tmpPosition[1], tmpPosition[2]);
            initialCameraViewUp_[renderer] = tmpPosition;

        }

        renderer->SetActiveCamera(camera);
        renderer->ResetCamera();
        displayRange_[renderer] = Range3D<indexType>();

        if (config_.renderBonds == true) {
            if (i == 3)
                renderBonds(renderer, config_.visualizeBondStrain, config_.refLatticeConstant);
            else
                renderBonds(renderer);
        }

        renderAtoms(renderer);

        if ((i == 3) && (config_.visualizeBondStrain == true))
            renderAbsStrainPlanes(renderer, *simbox_, strainField);

        auto le = largestExtent(bounds_);

        // autoscale individual views to fill the
        // individual render windows
        // camera 3 gets auto scale through ResetCamera function.
        if (i<3) camera->SetParallelScale(1.0*le);

        renderer->ResetCamera();
        
    }
}

//------------------------------------------------------------------------------

void SimulationBoxRenderer::addInterfaces(double commandLineHeightShare,
                                          const std::vector<Range3D<indexType>> interfaces)
{

    uint32_t j = 0;
    double yfrac = (1.0-commandLineHeightShare)/interfaces.size();
    std::vector<Vector3D<double>> cameraPositionsif;
    std::vector<double> xminsif, xmaxsif, yminsif, ymaxsif;
    Vector3D<spaceType> origin = {0.0, 0.0, 0.0};
    Vector3D<spaceType> spaceDiagonal = simbox_->getSize();
    
    auto largestExtent = [](double *boundsArray) {
        double result = 0.0;

        for (auto i: {0, 1, 2}) {
            if ((boundsArray[i*2+1] - boundsArray[i*2]) > result)
                result = (boundsArray[i*2+1] - boundsArray[i*2]);
        }

        return result;
    };

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    for (auto interface: interfaces) {
        xminsif.push_back(0.75);
        xmaxsif.push_back(1.0);

        yminsif.push_back(commandLineHeightShare+j*yfrac);
        ymaxsif.push_back(commandLineHeightShare+(j+1)*yfrac);

        Vector3D<double> tmpCameraPos = {};
        tmpCameraPos[0] = (spaceDiagonal[0] - origin[0]) / 2.0;
        tmpCameraPos[1] = (spaceDiagonal[1] - origin[1]) / 2.0;
        tmpCameraPos[2] = (spaceDiagonal[2] - origin[2]) * -5.0;

        cameraPositionsif.push_back(tmpCameraPos);

        j++;
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    j = 0;
    for (auto interface: interfaces) {
        Vector3D<spaceType> tmpPosition;
        vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();

        renderWindow_->AddRenderer(renderer);
        renderer->SetViewport(xminsif[j], yminsif[j], xmaxsif[j], ymaxsif[j]);
        renderer->SetBackground(1.0, 1.0, 1.0);

        vtkSmartPointer<vtkCamera> camera = vtkSmartPointer<vtkCamera>::New();
        tmpPosition  = {cameraPositionsif[j][0], cameraPositionsif[j][1], cameraPositionsif[j][2]};
        camera->SetPosition(tmpPosition[0], tmpPosition[1], tmpPosition[2]);
        initialCameraPosition_[renderer] = tmpPosition;
        
        vtkSmartPointer<vtkLegendBoxActor> legend = vtkSmartPointer<vtkLegendBoxActor>::New();
        vtkSmartPointer<vtkNamedColors> colors = vtkSmartPointer<vtkNamedColors>::New();
        double color[4];

        camera->ParallelProjectionOn();
        //renderer->InteractiveOff();

        legend->SetNumberOfEntries(1);

        legend->GetPositionCoordinate()->SetCoordinateSystemToView();
        legend->GetPositionCoordinate()->SetValue(0.5, -1.0);

        legend->GetPosition2Coordinate()->SetCoordinateSystemToView();
        legend->GetPosition2Coordinate()->SetValue(1.0, -0.5);

        legend->UseBackgroundOff();
        legend->BorderOff();

        const Atom* atom= simbox_->getLattice().getFirstAtomInLayer(interface.start[2], 2);
        camera->SetFocalPoint(cameraPositionsif[j][0], cameraPositionsif[j][1], atom->getPosition()[2]);
        camera->SetParallelScale(10.0);

        if ((j%2) == 0)
            renderer->SetBackground(1.0, 1.0, 1.0);
        else
            renderer->SetBackground(0.8, 0.8, 0.8);

        tmpPosition = {-1.0, 0.0, 0.0};
        camera->SetViewUp(tmpPosition[0], tmpPosition[1], tmpPosition[2]);
        initialCameraViewUp_[renderer] = tmpPosition;
        
        std::stringstream label{};
        label << "IF " << interface.start[2] << " (0 0 1)";

        colors->GetColor("black", color);
        legend->SetEntry(0, (vtkPolyData *)NULL, label.str().c_str(), color);

        renderer->AddActor(legend);

        renderer->SetActiveCamera(camera);
        renderer->ResetCamera();
        displayRange_[renderer] = interface;
        
        if (config_.renderBonds == true) {
            renderBonds(renderer);
        }

        renderAtoms(renderer);
        
        auto le = largestExtent(bounds_);

        // autoscale individual views to fill the
        // individual render windows
        // camera 3 gets auto scale through ResetCamera function.
        camera->SetParallelScale(1.0*le);

        renderer->ResetCamera();
        j++;

    }
}

//------------------------------------------------------------------------------

void SimulationBoxRenderer::renderPerspectiveWithInterfaces (SimulationBox &simbox,
                                                             const std::vector<Range3D<indexType>> interfaces,
                                                             const MaterialCollection &materials,
                                                             const Configuration &config,
                                                             std::shared_ptr<Field3D<double>> strainField,
                                                             bool persist) {
    config_ = config;
    simbox_ = &simbox;
    materials_ = &materials;

    double commandLineHeightShare = 0.04;

    if (strainField == nullptr){
        StrainTool st(*simbox_);
        strainField = st.getStrainField(*materials_, -1.0);
    }
        
    CLOG(TRACE, logName_) << "starting renderPerspectiveWithInterfaces()";

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    CLOG(TRACE, logName_) << "generating views of complete crystal and all interfaces";

    addCrystalViews(commandLineHeightShare, strainField);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    CLOG(TRACE, logName_) << "generating views of each single interface";
    
    addInterfaces(commandLineHeightShare, interfaces);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    CLOG(TRACE, logName_) << "generating command line";
    
    // add command line at bottom of window
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    
    renderWindow_->AddRenderer(renderer);
    renderer->SetViewport(0, 0, 1.0, commandLineHeightShare);
    renderer->SetBackground(1.0, 1.0, 1.0);
    
    vtkNamedColors *colors = vtkNamedColors::New();
    double color[3];
    colors->GetColor("black", color);
    
    vtkSmartPointer<vtkTextActor> textActor = vtkSmartPointer<vtkTextActor>::New();
        
    textActor->SetDisplayPosition(10, 4);

    textActor->GetTextProperty()->SetColor(color);
    textActor->GetTextProperty()->SetFontSize ( 20 );

    renderer->AddActor(textActor);
    KeyPressInteractorStyle::SafeDownCast(renderWindowInteractor_->
                                          GetInteractorStyle())->setCommandLine(textActor);

    renderWindow_->SetSize(windowWidth_,windowHeight_);
    renderWindow_->Render();

    if (persist == true) {
        CLOG(TRACE, logName_) << "starting Interactor";
        renderWindowInteractor_->Start();
    }
}

//******************************************************************************

void SimulationBoxRenderer::saveImage (const std::string filename) const {
    vtkSmartPointer<vtkWindowToImageFilter> w2iFilter = vtkSmartPointer<vtkWindowToImageFilter>::New();
    w2iFilter->SetInput(renderWindow_);

    std::string extension = filename.substr(filename.find_last_of(".")+1);
    std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);

    if ((extension == "jpg") || (extension == "jpeg")) {
        vtkSmartPointer<vtkJPEGWriter> jpgWriter = vtkSmartPointer<vtkJPEGWriter>::New();
        jpgWriter->SetInputConnection(w2iFilter->GetOutputPort());
        jpgWriter->SetFileName((config_.outputPreamble+filename).c_str());
        jpgWriter->Write();
    } else if (extension == "png") {
        vtkSmartPointer<vtkPNGWriter> pngWriter = vtkSmartPointer<vtkPNGWriter>::New();
        pngWriter->SetInputConnection(w2iFilter->GetOutputPort());
        pngWriter->SetFileName((config_.outputPreamble+filename).c_str());
        pngWriter->Write();
    } else {
        // \todo Raise exception
    }
}

//------------------------------------------------------------------------------

//! Determines, which layers (exclusively) contain modified atoms
//! The obtained information is stored in the member variables
//! modifiedAtomsInLayer and modifiedAtomsOnlyInLayer
//! \param modificationIndexLimit defines a maximum modification index to be considered
void SimulationBoxRenderer::evaluateModifiedAtomsInLayers (uint32_t modificationIndexLimit) {
    uint8_t outOfPlaneDimension = simbox_->getOutOfPlaneDimension();
    Vector3D<indexType> latticeSize = simbox_->getLattice().getSize();
    
    // reset default values first
    modifiedAtomsInLayer_.clear();
    modifiedAtomsInLayer_.resize(latticeSize[outOfPlaneDimension], false);
    modifiedAtomsOnlyInLayer_.clear();
    modifiedAtomsOnlyInLayer_.resize(latticeSize[outOfPlaneDimension],false);


    for (indexType i = 0; i < latticeSize[outOfPlaneDimension]; i++) {
        bool tmpLayerFullyModified = true;

        for (auto atom: simbox_->getLattice().getAtomsInLayer(i,outOfPlaneDimension)) {
            if (config_.filterModificationState.size() > 0) {
                if (atom->getState().any() == true) {
                    if (modificationIndexLimit > 0){
                        if ( ! (atom->getModificationOrder() < modificationIndexLimit) ) 
                            continue;
                    }

                    // if any atom was modified --> set dirty flag
                    modifiedAtomsInLayer_[i] = true;                
                } else {
                    tmpLayerFullyModified = false;
                }
            } else {
                if (atom->wasModified() == true){
                    if (modificationIndexLimit > 0){
                        if ( ! (atom->getModificationOrder() < modificationIndexLimit) )
                            continue;
                    }
                    
                    // if any atom was modified --> set dirty flag
                    modifiedAtomsInLayer_[i] = true;                
                        
                } else
                    tmpLayerFullyModified = false;
            }

        }

        modifiedAtomsOnlyInLayer_[i] = tmpLayerFullyModified;
    }
}

//------------------------------------------------------------------------------

//! Checks, if modified flags were evaluated already
//! \return true, if information on modified atoms in layers was set already
bool SimulationBoxRenderer::modifiedAtomsSet () const {
    bool result = false;

    if ((modifiedAtomsInLayer_.size() > 0) && (modifiedAtomsOnlyInLayer_.size() > 0))
        result = true;

    return result;
}

//------------------------------------------------------------------------------

uint32_t SimulationBoxRenderer::getMaxModificationIndex()
{
    return simbox_->getLattice().getMaxModificationIndex();
}

//------------------------------------------------------------------------------

void SimulationBoxRenderer::renderAtomsWithinLimit(vtkRenderer *renderer, const int limit)
{
    bool modifiedSave = config_.modifiedOnly;
    config_.modificationIndexMax = limit;
    config_.modifiedOnly = true;
    
    renderAtoms(renderer);
    renderBonds(renderer);
    
    config_.modifiedOnly = modifiedSave;    
}

//------------------------------------------------------------------------------

void SimulationBoxRenderer::resetCamera(vtkRenderer *renderer)
{

    vtkCamera *camera = renderer->GetActiveCamera();
    const Vector3D<spaceType> &position = initialCameraPosition_[renderer];
    camera->SetPosition(position[0], position[1], position[2]);

    const Vector3D<spaceType> &viewUp = initialCameraViewUp_[renderer];
    camera->SetViewUp(viewUp[0], viewUp[1], viewUp[2]);

    renderer->ResetCamera();
}

//##############################################################################

// Callback for the interaction
// class vtkMyCallback : public vtkCommand
// {
// public:
//   static vtkMyCallback *New()
//     { return new vtkMyCallback; }
//   void Execute(vtkObject *caller, unsigned long, void*) VTK_OVERRIDE
//   {
//       vtkRenderer *renderer = reinterpret_cast<vtkRenderer*>(caller);
//       vtkRenderWindowInteractor *iren = reinterpret_cast<vtkRenderWindowInteractor*>(caller);
//       auto key = iren->GetKeySym();
//       cout << "key " << key << " was pressed\n";
//   }
// };

vtkStandardNewMacro(KeyPressInteractorStyle);

//------------------------------------------------------------------------------

void KeyPressInteractorStyle::OnLeftButtonDown()
{
    vtkRenderWindowInteractor *rwi = this->Interactor;
    double clickBorders[4];

    CLOG(DEBUG, logName_) << "left button pressed at (" << rwi->GetEventPosition()[0]
         << "/" << rwi->GetEventPosition()[1] << ")";

    vtkRenderer *clickRenderer = rwi->FindPokedRenderer(rwi->GetEventPosition()[0],
                                                        rwi->GetEventPosition()[1]);
    clickRenderer->GetViewport(clickBorders);
    
    // only if clicked on command line at bottom of window
    if (clickBorders[1] == 0.0){
        commandLineMode_ = true;
        updateCommandLine();
    }

    // Forward events    
    vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
}

//------------------------------------------------------------------------------

void KeyPressInteractorStyle::OnChar()
{
    // Get the keypress
    vtkRenderWindowInteractor *rwi = this->Interactor;
    //    std::string key = rwi->GetKeySym();
    char key = rwi->GetKeyCode();

    CLOG(DEBUG, logName_) << "got character " << key;

    if (commandLineMode_){
        //        if (! (std::isalpha(key) || (key==' ') ) ) return;
        if (! std::isprint(key) ) return;
    
        commandString_.push_back(rwi->GetKeyCode());
        updateCommandLine();
        return;
    }
    
    if (key == 's'){
        // switch renderers to middle window

        double clickBorders[4], midBorders[4];
        double clickBackground[3], midBackground[3];

        vtkRenderer *clickRenderer = rwi->FindPokedRenderer(rwi->GetEventPosition()[0],
                                                            rwi->GetEventPosition()[1]);
        vtkRenderer *midRenderer = rwi->FindPokedRenderer(windowWidth_ /2, windowHeight_ /2);

        clickRenderer->GetViewport(clickBorders);
        midRenderer->GetViewport(midBorders);

        // Do not switch areas on the left border
        if (clickBorders[0] != 0.0){

            clickRenderer->GetBackground(clickBackground);
            midRenderer->GetBackground(midBackground);

            midRenderer->SetViewport(clickBorders);
            clickRenderer->SetViewport(midBorders);

            clickRenderer->SetBackground(midBackground);
            midRenderer->SetBackground(clickBackground);
        }

        rwi->Render();
    }
    else if (key == 'c'){

        vtkRenderer *clickRenderer = rwi->FindPokedRenderer(rwi->GetEventPosition()[0],
                                                            rwi->GetEventPosition()[1]);

        simboxRenderer_->resetCamera(clickRenderer);
        rwi->Render();
    }
    else if (key == 'u'){

        commandLineMode_ = true;
        updateCommandLine ();
    }
    else if (key == 'b'){
        toggleBonds(rwi);
    }
    else if (key == 'a'){
        animate(rwi);

    }
    
    // Forward events
    vtkInteractorStyleTrackballCamera::OnChar();
}

//------------------------------------------------------------------------------

void KeyPressInteractorStyle::runCommand(const std::string command)
{
    CLOG(DEBUG, logName_) << "Running command \'" << command << "\'";
    std::string commandRef;
    
    if (command == "bonds"){
        toggleBonds(this->Interactor);
        return;
    }

    if (command == "animate"){
        animate(this->Interactor);
        return;
    }

    commandRef = "timeout ";
    if ( (command.length() > commandRef.length()) &&
         (command.substr(0,commandRef.length()) == commandRef) ){
        animationTimeout_ = std::stoi(command.substr(commandRef.length()));
        return;
    }

    commandRef = "image-preamble ";    
    if ( (command.length() > commandRef.length() ) &&
         (command.substr(0,commandRef.length() ) == commandRef) ){
        imagePreamble_ = command.substr(commandRef.length());
        return;
    }
    //this->Interactor->GetRenderWindow()->MakeCurrent();
}

//------------------------------------------------------------------------------

void KeyPressInteractorStyle::OnKeyPress()
{
    // Get the keypress
    vtkRenderWindowInteractor *rwi = this->Interactor;
    std::string key = rwi->GetKeySym();

    CLOG(DEBUG, logName_) << "key " << key  << " was pressed";
    
    if (key == "Return"){
        if (commandLineMode_ == true) {
            runCommand(commandString_);
            commandString_.clear();
            commandLineMode_ = false;
        }
        else
            commandLineMode_ = true;
        
        updateCommandLine();
    }
    if (key == "BackSpace"){
        if (! commandString_.empty()) {
            commandString_.pop_back();
            //            cout << "after removal " << commandString_ << std::endl;
        }
        updateCommandLine();
    }
    if (key == "Escape"){
        if (commandLineMode_){
            commandString_.clear();
            commandLineMode_ = false;
            updateCommandLine();
        }
        else
            rwi->TerminateApp();
    }

    // Forward events
    vtkInteractorStyleTrackballCamera::OnKeyPress();
}

//------------------------------------------------------------------------------

void KeyPressInteractorStyle::setSize(const uint16_t windowWidth, const uint16_t windowHeight)
{
    windowWidth_ = windowWidth;
    windowHeight_ = windowHeight;
}

//------------------------------------------------------------------------------

void KeyPressInteractorStyle::setSimboxRenderer(SimulationBoxRenderer *simboxRenderer)
{
    simboxRenderer_ = simboxRenderer;
}

//------------------------------------------------------------------------------

void KeyPressInteractorStyle::setLogName(std::string logName)
{
    logName = logName_;
}

//------------------------------------------------------------------------------

void KeyPressInteractorStyle::setCommandLine(vtkTextActor *commandLineActor)
{
    commandLineActor_ = commandLineActor;
    updateCommandLine();
}

//------------------------------------------------------------------------------

void KeyPressInteractorStyle::updateCommandLine()
{

    std::string displayString = "Command $ " + commandString_;

    if (commandLineMode_)
        displayString += '_';

    commandLineActor_->SetInput(displayString.c_str());
    this->Interactor->Render();
}

//------------------------------------------------------------------------------

//! \todo find better property to identify actors that shall be deleted
void KeyPressInteractorStyle::toggleBonds(vtkRenderWindowInteractor *rwi)
{

    vtkRenderer *ren;
    vtkRendererCollection *rc = rwi->GetRenderWindow()->GetRenderers();
    vtkActor *actor;
    vtkCollectionSimpleIterator rsit;
    int actorCount;

    for (rc->InitTraversal(rsit); (ren = rc->GetNextRenderer(rsit)); ){

        vtkActorCollection *ac = ren->GetActors();
        vtkCollectionSimpleIterator ait;
        actorCount = 0;
        
        for ( ac->InitTraversal(ait); (actor=ac->GetNextActor(ait)); ){
            vtkProperty *prop = actor->GetProperty();

            // currently actor is identified by having line width == 2
            if (prop->GetLineWidth() == 2){
                ren->RemoveActor(actor);
                actorCount ++;
            }
        }

        if (actorCount == 0)
            simboxRenderer_->renderBonds(ren);
    }
    rwi->Render();
}

//------------------------------------------------------------------------------

void KeyPressInteractorStyle::animate(vtkRenderWindowInteractor *rwi)
{

    vtkRenderer *ren;
    vtkRendererCollection *rc = rwi->GetRenderWindow()->GetRenderers();
    vtkActor *actor;
    vtkCollectionSimpleIterator rsit;
    int actorCount;

    uint32_t maxModificationIndex = simboxRenderer_->getMaxModificationIndex();

    for (int i=1; i < maxModificationIndex; i++){
        
        for (rc->InitTraversal(rsit); (ren = rc->GetNextRenderer(rsit)); ){

            vtkActorCollection *ac = ren->GetActors();
            vtkCollectionSimpleIterator ait;
        
            for ( ac->InitTraversal(ait); (actor=ac->GetNextActor(ait)); ){
                vtkProperty *prop = actor->GetProperty();

                ren->RemoveActor(actor);

            }
            simboxRenderer_->renderAtomsWithinLimit(ren, i);

        }

        rwi->Render();
        std::stringstream name;
        name << imagePreamble_ << "animate_" << std::right <<
            std::setw(3) << std::setfill('0') << i << ".png";
        simboxRenderer_->saveImage(name.str());
        std::this_thread::sleep_for(std::chrono::milliseconds(animationTimeout_));
    }
}


//------------------------------------------------------------------------------

// void KeyPressInteractorStyle::OnMouseMove()
// {
//     int x, y;

//     this->Interactor->GetMousePosition(&x, &y);
//     if (x > (1280/4)) vtkInteractorStyleTrackballCamera::OnMouseMove();
// }

//##############################################################################

#endif

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
