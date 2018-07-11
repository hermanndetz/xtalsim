/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/

#ifndef __XY_PLOT_H__
#define __XY_PLOT_H__

#include "projectConfigure.h"

#ifdef __VTK__

#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkTable.h>

#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkTable.h>
#include <vtkVariantArray.h>
#include <vtkXYPlotActor.h>
#include <vtkChartXY.h>
#include <vtkPlot.h>
#include <vtkContextView.h>
#include <vtkContextScene.h>
#include <vtkPen.h>
#include <vtkAxis.h>

#include <misc/CsvHandler.h>
#include <visualization/DataSeries.h>

typedef std::tuple<double,double> axesLimits;

class XYPlot {
    private:
    const std::string logName_;

        std::vector<vtkSmartPointer<vtkTable>> tables_;

        std::string xLabel_="x Axis";
        std::string yLabel_="y Axis";
        std::string title_;

        axesLimits xLimits_;
        axesLimits yLimits_;
        bool xLimitsSet=false;
        bool yLimitsSet=false;
        bool xLimitFixed=false;
        bool yLimitFixed=false;

        vtkSmartPointer<vtkRenderer> renderer_;
        vtkSmartPointer<vtkRenderWindow> renderWindow_;
        vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor_;

    public:
        //! Constructor
        XYPlot (std::string title="", const char *logName="XYPlot");
        XYPlot (vtkSmartPointer<vtkTable> _table, std::string title="", const char *logName="XYPlot");
        XYPlot (std::vector<vtkSmartPointer<vtkTable>> tables, std::string title="", const char *logName="XYPlot");

        //! Destructor
        ~XYPlot ();

        void    plot ();
        void    plot (std::vector<DataSeries> dataSeries);
        void    plot (const std::string fileName);
        void    setTable (vtkSmartPointer<vtkTable> _table);
        void    setAxesLimits (axesLimits &xLimits, axesLimits &yLimits, bool xFixed=true, bool yFixed=true);
        void    setAxesLabels (std::string xAxisLabel="x Axis", std::string yAxisLabel="y Axis");
};

#endif
#endif


// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
