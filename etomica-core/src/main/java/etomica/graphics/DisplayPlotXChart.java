package etomica.graphics;

import etomica.data.DataSet;
import etomica.data.DataSetListener;
import etomica.data.IDataInfo;
import etomica.data.IDataSink;
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
    private final List<XYSeries> series;
    private final Map<String, Unit> unitMap;
    private Unit defaultUnit = null;
    private Unit xUnit = null;

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
        this.series.add(series);
        return sink;
    }

    public XYSeries getSeries(String name) {
        return this.plot.getSeriesMap().get(name);
    }

    public void setUnit(String seriesName, Unit unit) {
        this.unitMap.put(seriesName, unit);
    }

    public void setDefaultUnit(Unit unit) {
        this.defaultUnit = unit;
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

    public XYChart getPlot() {
        return this.plot;
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

        Set<String> newNames = IntStream.range(0, dataSet.getDataCount())
                .mapToObj(this.dataSet::getName)
                .collect(Collectors.toSet());
        for (String s : this.plot.getSeriesMap().keySet()) {
            if (!newNames.contains(s)) {
                XYSeries ser = this.plot.getSeriesMap().remove(s);
                this.series.remove(ser);
                this.unitMap.remove(s);
            }
        }
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
                double[] xValues = dataInfo.getXDataSource().getIndependentData(0).getData();
                xValues = Arrays.stream(xValues).map(x -> this.xUnit.fromSim(x)).toArray();
                double[] data = ((DataFunction)dataSet.getData(i)).getData();
                data = Arrays.stream(data).map(yUnit::fromSim).toArray();

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
