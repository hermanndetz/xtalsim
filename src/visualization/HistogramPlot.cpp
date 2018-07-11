/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/

#ifdef __VTK__

#include "HistogramPlot.h"

//------------------------------------------------------------------------------

//! Constructor
HistogramPlot::HistogramPlot (const char *logName):
    logName_(logName)
{
}

//------------------------------------------------------------------------------

HistogramPlot::HistogramPlot (vtkSmartPointer<vtkTable> table, const char *logName):
    logName_(logName), table_(table)
{
}

//------------------------------------------------------------------------------

HistogramPlot::~HistogramPlot () {

}

//------------------------------------------------------------------------------

void HistogramPlot::plot (void) {
    if (renderer_ == nullptr)
        renderer_ = vtkSmartPointer<vtkRenderer>::New();

    if (renderWindow_ == nullptr)
        renderWindow_ = vtkSmartPointer<vtkRenderWindow>::New();

    if (renderWindowInteractor_ == nullptr)
        renderWindowInteractor_ = vtkSmartPointer<vtkRenderWindowInteractor>::New();

    vtkSmartPointer<vtkContextView> view = vtkSmartPointer<vtkContextView>::New();
    view->GetRenderer()->SetBackground(1.0, 1.0, 1.0);

    vtkSmartPointer<vtkChartXY> chart = vtkSmartPointer<vtkChartXY>::New();
    view->GetScene()->AddItem(chart);

    uint32_t i = 0;

    for (i = 1; i < table_->GetNumberOfColumns(); i++) {
        Color col = Color("FF0000");
        vtkPlot *line = chart->AddPlot(vtkChart::BAR);

        line->SetInputData(table_, 0, i);
        line->SetColor(col.get(Color::Red), col.get(Color::Green),
		       col.get(Color::Blue), 255);
        line->SetWidth(5);
        line->SetUseIndexForXSeries(false);
        line = chart->AddPlot(vtkChart::BAR);
    }

    /*
    for (auto item: dataSeries) {
        Color col = item.GetColor();
        vtkPlot *line = chart->AddPlot(vtkChart::LINE);
        //line->SetInputData(dTable, item.GetXIndex(), item.GetYIndex());
        line->SetInputData(tables_[item.GetTableIndex()], item.GetXIndex(), item.GetYIndex());
        line->SetColor(col.get(Color::Red), col.get(Color::Green),
		       col.get(Color::Blue), 255);
        line->SetWidth(1.5);
        line->SetUseIndexForXSeries(false);
        line = chart->AddPlot(vtkChart::LINE);
    }
    */


    view->GetInteractor()->Initialize();
    view->GetInteractor()->Start();

}

//------------------------------------------------------------------------------

//! Can be used to set the data table to be used for plots.
//! \param _table vtkSmartPointer to vtkTable.
void HistogramPlot::setTable (vtkSmartPointer<vtkTable> table) {
    table_ = table;
}

#endif
//##############################################################################


// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
