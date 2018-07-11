/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/

#include "XYPlot.h"

#ifdef __VTK__

//! Constructor
XYPlot::XYPlot (std::string title, const char *logName):
    title_(title), logName_(logName)
{
}

//------------------------------------------------------------------------------

//! Constructor
XYPlot::XYPlot (vtkSmartPointer<vtkTable> table, std::string title, const char *logName):
    title_(title), logName_(logName)
{
    tables_.push_back(table);
}

//------------------------------------------------------------------------------

//! Constructor
XYPlot::XYPlot (std::vector<vtkSmartPointer<vtkTable>> tables, std::string title, const char *logName):
    title_(title), logName_(logName)
{
    for (auto item: tables)
        tables_.push_back(item);
}

//------------------------------------------------------------------------------

//! Destructor
XYPlot::~XYPlot () {

}

//------------------------------------------------------------------------------

//! Plots the contents of table 0 (col 1 vs col 0).
void XYPlot::plot () {
    std::vector<DataSeries> dataSeries;

    Color col(0, 0, 0);
    DataSeries ds = DataSeries(0, 0, 1, col);
    dataSeries.push_back(ds);

    plot(dataSeries);
}

//------------------------------------------------------------------------------

//! Plots data series from a table, which is loaded.
//! \param dataSeries should contain a vector of DataSeries objects, which define the individual data series.
//!
//! \note This function still contains a few fragments, which use an alternative
//!       way to plot a function using VTK.
void XYPlot::plot (std::vector<DataSeries> dataSeries) {
    /*
    vtkSmartPointer<vtkXYPlotActor> plot = vtkSmartPointer<vtkXYPlotActor>::New();
    plot->ExchangeAxesOff();
    plot->SetLabelFormat("%g");
    plot->SetXTitle(xLabel_.c_str());
    plot->SetYTitle(yLabel_.c_str());
    plot->SetXValuesToIndex();

    plot->SetPlotColor(0,1,0,0);
    plot->SetPlotColor(1,0,1,0);
    */

    if (renderer_ == nullptr)
        renderer_ = vtkSmartPointer<vtkRenderer>::New();

    if (renderWindow_ == nullptr)
        renderWindow_ = vtkSmartPointer<vtkRenderWindow>::New();

    if (renderWindowInteractor_ == nullptr)
        renderWindowInteractor_ = vtkSmartPointer<vtkRenderWindowInteractor>::New();

    /* for plotting as field data
    renderer_->AddActor(plot);
    renderWindow_->AddRenderer(renderer);
    renderWindowInteractor_->SetRenderWindow(renderWindow);

    renderWindowInteractor_->Initialize();
    renderWindowInteractor_->Start();
    */


    vtkSmartPointer<vtkContextView> view = vtkSmartPointer<vtkContextView>::New();
    view->GetRenderer()->SetBackground(1.0, 1.0, 1.0);

    vtkSmartPointer<vtkChartXY> chart = vtkSmartPointer<vtkChartXY>::New();
    view->GetScene()->AddItem(chart);

    if (title_ != "")
        chart->SetTitle(title_);

    chart->GetAxis(1)->SetTitle(xLabel_.c_str());
    chart->GetAxis(0)->SetTitle(yLabel_.c_str());

    chart->SetShowLegend(true);
    
    if (xLimitsSet == true) {
        // x-axis has index 1 in VTK.
        chart->GetAxis(1)->SetRange(std::get<0>(xLimits_), std::get<1>(xLimits_));
        chart->GetAxis(1)->SetMinimumLimit(std::get<0>(xLimits_));
        chart->GetAxis(1)->SetMaximumLimit(std::get<1>(xLimits_));
        chart->GetAxis(1)->SetUnscaledMinimumLimit(std::get<0>(xLimits_));
        chart->GetAxis(1)->SetUnscaledMaximumLimit(std::get<1>(xLimits_));

        if (xLimitFixed == true)
            chart->GetAxis(1)->SetBehavior(vtkAxis::FIXED);
    }
    
    if (yLimitsSet == true) {
        // y-axis has index 0 in VTK.
        chart->GetAxis(0)->SetRange(std::get<0>(yLimits_), std::get<1>(yLimits_));
        chart->GetAxis(0)->SetMinimumLimit(std::get<0>(yLimits_));
        chart->GetAxis(0)->SetMaximumLimit(std::get<1>(yLimits_));

        if (yLimitFixed == true)
            chart->GetAxis(0)->SetBehavior(vtkAxis::FIXED);
    }
    
    for (auto item: dataSeries) {
        Color col = item.GetColor();
        vtkPlot *line = chart->AddPlot(vtkChart::LINE);
        //line->SetInputData(dTable, item.GetXIndex(), item.GetYIndex());

        line->SetInputData(tables_[item.GetTableIndex()], item.GetXIndex(), item.GetYIndex());
        line->SetColor(col.get(Color::Red), col.get(Color::Green),
		       col.get(Color::Blue), 255);
        line->SetWidth(2.5);
        line->SetUseIndexForXSeries(false);
    }

    /*
    line->SetInputData(dTable, 0, 2);
    line->SetColor(255, 0, 0, 255);
    line->SetWidth(2.0);
    line->GetPen()->SetLineType(vtkPen::DASH_LINE);
    */

    view->GetInteractor()->Initialize();
    view->GetInteractor()->Start();
}

//------------------------------------------------------------------------------

//! Plots the first data series (col2 vs col1) of a given .csv or .tsv file.
//! \param fileName Path to the input file.
void XYPlot::plot (const std::string fileName) {
    tables_.push_back(CsvHandler(fileName, "CSVFile").get());

    plot();
}

//------------------------------------------------------------------------------

//! Can be used to set the data table to be used for plots.
//! \param _table vtkSmartPointer to vtkTable.
void XYPlot::setTable (vtkSmartPointer<vtkTable> table) {
    tables_.push_back(table);
}

//------------------------------------------------------------------------------

//! Can be used to set axes limits.
//! \param xLimit tuple defining limits for x axis.
//! \param yLimit tuple defining limits for y axis.
//! \param xFixed disables scaling of x axis if true (default)
//! \param yFixed disables scaling of y axis if true (default)
//! If called, this function enforces the limits per default.
//! \todo This basic implementation always sets limits for both axes.
void XYPlot::setAxesLimits (axesLimits &xLimits, axesLimits &yLimits, 
        bool xFixed, bool yFixed) {
    xLimits_ = xLimits;
    yLimits_ = yLimits;
    xLimitFixed = xLimitFixed;
    yLimitFixed = yLimitFixed;
    xLimitsSet = true;
    yLimitsSet = true;
}

//------------------------------------------------------------------------------

void XYPlot::setAxesLabels (std::string xAxisLabel, std::string yAxisLabel) {

    xLabel_ = xAxisLabel;
    yLabel_ = yAxisLabel;
}

#endif
//##############################################################################


// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
