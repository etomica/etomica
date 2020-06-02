package etomica.graphics;

import etomica.data.DataSet;
import etomica.data.DataSetListener;
import org.knowm.xchart.XChartPanel;
import org.knowm.xchart.XYChart;
import org.knowm.xchart.XYChartBuilder;

import javax.swing.*;

public class DisplayPlotXChart extends Display implements DataSetListener {
    private final DataSet dataSet;
    private final XYChart plot;
    private final JPanel panel;

    public DisplayPlotXChart(DataSet dataSet) {
        this.dataSet = dataSet;
        this.plot = new XYChartBuilder()
                .build();
        this.panel = new JPanel();
        XChartPanel<XYChart> panel = new XChartPanel<>(this.plot);
        this.panel.add(panel);
    }

    @Override
    public void dataChanged(DataSet dataSet) {

    }

    @Override
    public void dataCountChanged(DataSet dataSet) {

    }

    private void update() {

    }
}
