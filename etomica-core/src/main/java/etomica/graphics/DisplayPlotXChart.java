package etomica.graphics;

import etomica.data.DataSet;
import etomica.data.DataSetListener;
import etomica.data.types.DataFunction;
import org.knowm.xchart.*;

import javax.swing.*;
import java.awt.*;
import java.util.ArrayList;
import java.util.List;

public class DisplayPlotXChart extends Display implements DataSetListener {
    private final DataSet dataSet;
    private final XYChart plot;
    private final XChartPanel<XYChart> panel;
    private final List<XYSeries> series;

    public DisplayPlotXChart() {
        this(new DataSet());
    }

    public DisplayPlotXChart(DataSet dataSet) {
        this.dataSet = dataSet;
        this.dataSet.addDataListener(this);
        this.plot = new XYChartBuilder()
                .build();
        this.panel = new XChartPanel<>(this.plot);
        this.series = new ArrayList<>();
    }

    public Component graphic() {
        return this.panel;
    }

    public void repaint() {
        panel.repaint();
    }

    public DataSet getDataSet() {
        return this.dataSet;
    }

    @Override
    public void dataChanged(DataSet dataSet) {
        this.update();
    }

    @Override
    public void dataCountChanged(DataSet dataSet) {
        for (String s : this.plot.getSeriesMap().keySet()) {
            this.plot.removeSeries(s);
        }
        this.series.clear();

        for (int i = 0; i < this.dataSet.getDataCount(); i++) {
            if (this.dataSet.getDataInfo(i) instanceof DataFunction.DataInfoFunction) {
                double[] xValues = ((DataFunction.DataInfoFunction) dataSet.getDataInfo(i)).getXDataSource().getIndependentData(0).getData();
                double[] data = ((DataFunction)dataSet.getData(i)).getData();

                XYSeries series = this.plot.addSeries("Series" + i, xValues, data);
                this.series.add(series);
            }
        }
    }

    private void update() {
        int nSeries = this.dataSet.getDataCount();

        for (int i = 0; i < nSeries; i++) {
            if (this.dataSet.getDataInfo(i) instanceof DataFunction.DataInfoFunction) {
                double[] xValues = ((DataFunction.DataInfoFunction) dataSet.getDataInfo(i)).getXDataSource().getIndependentData(0).getData();
                double[] data = ((DataFunction)dataSet.getData(i)).getData();

                XYSeries series = this.series.get(i);
                this.plot.updateXYSeries(series.getName(), xValues, data, null);
            }
        }
        if (this.panel.isShowing()) {
            this.panel.revalidate();
            this.panel.repaint();
        }
    }
}
