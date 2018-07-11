/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/

#ifndef __HISTOGRAM_PLOT_H__
#define __HISTOGRAM_PLOT_H__

#include "projectConfigure.h"

#ifdef __VTK__

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkTable.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkContextView.h>
#include <vtkContextScene.h>
#include <vtkChartXY.h>
#include <vtkPlot.h>

#include <misc/Color.h>

class HistogramPlot {
    private:
    const std::string logName_;

        vtkSmartPointer<vtkTable> table_;

        vtkSmartPointer<vtkRenderer> renderer_;
        vtkSmartPointer<vtkRenderWindow> renderWindow_;
        vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor_;

    public:
        //! Constructor
        HistogramPlot (const char *logName="HistogramPlot");
        HistogramPlot (vtkSmartPointer<vtkTable> _table, const char *logName="HistogramPlot");

        //! Destructor
        ~HistogramPlot ();

        void plot (void);
        void setTable (vtkSmartPointer<vtkTable> table);
};

#endif

#endif


// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
