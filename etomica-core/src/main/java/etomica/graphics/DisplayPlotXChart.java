package etomica.graphics;

import etomica.data.*;
import etomica.data.types.DataFunction;
import etomica.units.Unit;
import etomica.units.dimensions.Null;
import etomica.units.systems.UnitSystem;
import org.knowm.xchart.*;
import org.knowm.xchart.internal.chartpart.SelectionZoom;

import javax.swing.*;
import java.awt.*;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class DisplayPlotXChart extends Display implements DataSetListener {
    private final DataSet dataSet;
    private final XYChart plot;
    private final XChartPanel<XYChart> panel;
    private final Map<String, Unit> unitMap;
    private Unit defaultUnit = null;
    private Unit xUnit = null;

    private final List<DataTagBag> labelList = new ArrayList<>();
    private final List<DataTagBag> unitList = new ArrayList<>();
    private final List<DataTagBag> drawLineList = new ArrayList<>();

    public DisplayPlotXChart() {
        this(new DataSet());
    }

    public DisplayPlotXChart(DataSet dataSet) {
        this.dataSet = dataSet;
        this.dataSet.addDataListener(this);
        this.plot = new XYChartBuilder()
                .build();
        this.panel = new XChartPanel<>(this.plot);
        this.unitMap = new HashMap<>();
    }

    public Component graphic() {
        return this.panel;
    }

    public DataSet getDataSet() {
        return this.dataSet;
    }

    public IDataSink makeSink(String name) {
        IDataSink sink = this.dataSet.makeDataSink(name);
        XYSeries series = this.plot.addSeries(name, new double[]{Double.NaN}, new double[]{Double.NaN});
        this.plot.updateXYSeries(name, new double[0], new double[0], null);
        return sink;
    }

    public void removeSeries(String name) {
        this.plot.removeSeries(name);
        this.unitMap.remove(name);
    }

    public XYSeries getSeries(String name) {
        return this.plot.getSeriesMap().get(name);
    }

    public void setUnit(String seriesName, Unit unit) {
        this.unitMap.put(seriesName, unit);
    }

    public void setUnit(DataTag[] dataTags, Unit newUnit) {
        this.unitList.add(new DataTagBag(dataTags, newUnit));
    }

    public void setDefaultUnit(Unit unit) {
        this.defaultUnit = unit;
    }

    public void setUnit(Unit unit) {
        this.setDefaultUnit(unit);
    }

    public void setXUnit(Unit unit) {
        this.xUnit = unit;
    }

    public void setXLabel(String label) {
        this.plot.setXAxisTitle(label);
    }

    public void setYLabel(String label) {
        this.plot.setYAxisTitle(label);
    }

    public void setXLog(boolean log) {
        this.plot.getStyler().setXAxisLogarithmic(log);
    }

    public void setYLog(boolean log) {
        this.plot.getStyler().setYAxisLogarithmic(log);
    }

    public void setDoLegend(boolean legend) {
        this.plot.getStyler().setLegendVisible(legend);

        for (int i = 0; i < this.dataSet.getDataCount(); i++) {
            String dataLabel = this.dataSet.getDataInfo(i).getLabel();
            DataTagBag tagLabel = DataTagBag.getDataTagBag(this.labelList, this.dataSet.getDataInfo(i).getTags());
            if (tagLabel != null) {
                dataLabel = (String)tagLabel.object;
            }
            this.getSeries(this.dataSet.getName(i)).setLabel(dataLabel);
        }
    }

    public void setDoDrawLines(DataTag[] dataTags, boolean doDrawLines) {
        this.drawLineList.add(new DataTagBag(dataTags, doDrawLines));
    }

    public void setLegend(DataTag[] dataTags, String label) {
        this.labelList.add(new DataTagBag(dataTags, label));
    }

    public void setTitle(String title) {
        this.plot.setTitle(title);
    }

    public XYChart getChart() {
        return this.plot;
    }

    public XChartPanel<XYChart> getPanel() {
        return this.panel;
    }

    public DisplayPlotXChart getPlot() {
        return this;
    }

    public Dimension getPreferredSize() {
        return this.panel.getPreferredSize();
    }

    @Override
    public void dataChanged(DataSet dataSet) {
        this.update();
    }

    @Override
    public void dataCountChanged(DataSet dataSet) {
        if (this.dataSet.getDataCount() > 0) {
            IDataInfo xDataInfo = ((DataFunction.DataInfoFunction) this.dataSet.getDataInfo(0)).getXDataSource().getIndependentDataInfo(0);
            if (this.plot.getXAxisTitle().isEmpty()) {
                this.plot.setXAxisTitle(xDataInfo.getLabel());
            }

            if (this.xUnit == null) {
                this.xUnit = xDataInfo.getDimension().getUnit(UnitSystem.SIM);
            }
        }
        this.setDoLegend(this.plot.getStyler().isLegendVisible());

        for (int i = 0; i < this.dataSet.getDataCount(); i++) {
            DataTagBag tagDrawLines = DataTagBag.getDataTagBag(this.drawLineList, this.dataSet.getDataInfo(i).getTags());
            if (tagDrawLines != null) {
                Boolean dataDrawLines = (Boolean)tagDrawLines.object;
                this.getSeries(this.dataSet.getName(i))
                        .setXYSeriesRenderStyle(dataDrawLines ? XYSeries.XYSeriesRenderStyle.Line : XYSeries.XYSeriesRenderStyle.Scatter);
            }
        }

        for (int i = 0; i < this.dataSet.getDataCount(); i++) {
            Unit dataUnit = this.defaultUnit;
            DataTagBag tagUnit = DataTagBag.getDataTagBag(this.unitList, this.dataSet.getDataInfo(i).getTags());
            if (tagUnit != null) {
                dataUnit = (Unit)tagUnit.object;
            }
            this.setUnit(this.dataSet.getName(i), dataUnit);
        }
    }

    public void doUpdate() {
        this.update();
    }

    private void update() {
        int nSeries = this.dataSet.getDataCount();

        for (int i = 0; i < nSeries; i++) {
            if (this.dataSet.getDataInfo(i) instanceof DataFunction.DataInfoFunction) {
                String seriesName = this.dataSet.getName(i);
                DataFunction.DataInfoFunction dataInfo = (DataFunction.DataInfoFunction) dataSet.getDataInfo(i);
                Unit yUnit = this.unitMap.computeIfAbsent(seriesName, name -> {
                    return this.defaultUnit == null ? dataInfo.getDimension().getUnit(UnitSystem.SIM) : defaultUnit;
                });
                double[] data = ((DataFunction)dataSet.getData(i)).getData();
                double[] xValues = dataInfo.getXDataSource().getIndependentData(0).getData();

                ArrayList<Double> filteredData = new ArrayList<>(data.length);
                ArrayList<Double> filteredXValues = new ArrayList<>(xValues.length);
                IntStream.range(0, xValues.length)
                        .filter(idx -> {
                            boolean xLogBad = this.plot.getStyler().isXAxisLogarithmic() && xValues[idx] <= 0;
                            boolean yLogBad = this.plot.getStyler().isYAxisLogarithmic() && data[idx] <= 0;
                            return !Double.isNaN(yUnit.fromSim(data[idx])) && !xLogBad && !yLogBad;
                        })
                        .forEach(idx -> {
                            filteredData.add(yUnit.fromSim(data[idx]));
                            filteredXValues.add(xUnit.fromSim(xValues[idx]));
                        });



                SwingUtilities.invokeLater(() -> {
                    this.plot.getSeriesMap().get(seriesName).setEnabled(!filteredData.isEmpty());
                    this.plot.updateXYSeries(seriesName, filteredXValues, filteredData, null);
                });
            }
        }
        if (this.panel.isShowing()) {
            this.panel.revalidate();
            this.panel.repaint();
        }
    }
}
