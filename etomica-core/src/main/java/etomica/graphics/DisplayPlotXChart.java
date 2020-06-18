package etomica.graphics;

import etomica.data.*;
import etomica.data.types.DataFunction;
import etomica.units.Unit;
import etomica.units.systems.UnitSystem;
import org.knowm.xchart.*;

import javax.swing.*;
import javax.swing.plaf.LayerUI;
import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.geom.Rectangle2D;
import java.util.*;
import java.util.List;

public class DisplayPlotXChart extends Display implements DataSetListener {
    private final DataSet dataSet;
    private final XYChart plot;
    private final XChartPanel<XYChart> panel;
    private final Map<String, Unit> unitMap;
    private Unit defaultUnit = null;
    private Unit xUnit = null;
    private final Map<String, double[]> seriesXValues = new HashMap<>();
    private final Map<String, double[]> seriesYValues = new HashMap<>();

    private final List<DataTagBag> labelList = new ArrayList<>();
    private final List<DataTagBag> unitList = new ArrayList<>();
    private final List<DataTagBag> drawLineList = new ArrayList<>();

    private final ZoomUI zoomUI = new ZoomUI();
    private final JLayer<XChartPanel<XYChart>> layer;

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
        ZoomMouseAdapter zoomAdapter = new ZoomMouseAdapter();
        this.panel.addMouseListener(zoomAdapter);
        this.panel.addMouseMotionListener(zoomAdapter);

        this.layer = new JLayer<>(this.panel, this.zoomUI);
    }

    public Component graphic() {
        return this.layer;
    }

    public DataSet getDataSet() {
        return this.dataSet;
    }

    public IDataSink makeSink(String name) {
        IDataSink sink = this.dataSet.makeDataSink(name);
        XYSeries series = this.plot.addSeries(name, new double[]{Double.NaN}, new double[]{Double.NaN});
        series.setEnabled(false);
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

    public void setXRange(double min, double max) {
        this.plot.getStyler().setXAxisMin(min).setXAxisMax(max);
    }

    public void setYRange(double min, double max) {
        this.plot.getStyler().setYAxisMin(min).setYAxisMax(max);
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

    public JComponent getPanel() {
        return this.layer;
    }

    public void setSize(int w, int h) {
        this.panel.setSize(w, h);
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

    private static double[] getSeriesArray(Map<String, double[]> map, String name, int size) {
        return map.compute(name, (k, arr) -> {
            if (arr == null || arr.length != size) {
                return new double[size];
            } else {
                return arr;
            }
        });

    }

    private void update() {
        int nSeries = this.dataSet.getDataCount();

        double[][] newXValues = new double[nSeries][];
        double[][] newYValues = new double[nSeries][];

        boolean xAxisLog = this.plot.getStyler().isXAxisLogarithmic();
        boolean yAxisLog = this.plot.getStyler().isYAxisLogarithmic();

        for (int i = 0; i < nSeries; i++) {
            if (this.dataSet.getDataInfo(i) instanceof DataFunction.DataInfoFunction) {
                String seriesName = this.dataSet.getName(i);
                DataFunction.DataInfoFunction dataInfo = (DataFunction.DataInfoFunction) dataSet.getDataInfo(i);

                Unit yUnit = this.unitMap.computeIfAbsent(
                        seriesName,
                        name -> this.defaultUnit == null ? dataInfo.getDimension().getUnit(UnitSystem.SIM) : defaultUnit);

                double[] data = ((DataFunction)dataSet.getData(i)).getData();
                double[] xValues = dataInfo.getXDataSource().getIndependentData(0).getData();

                int leadingNans = 0;
                boolean leadingNansDone = false;
                int trailingNans = 0;
                for (int j = 0; j < xValues.length; j++) {
                    double x = xValues[j];
                    double y = data[j];

                    boolean badValue = Double.isNaN(x) || Double.isNaN(y) || (xAxisLog && x <= 0) || (yAxisLog && y <= 0);
                    if (!leadingNansDone && badValue) {
                        leadingNans++;
                    } else {
                        leadingNansDone = true;
                    }

                    if (leadingNansDone && badValue) {
                        trailingNans++;
                    } else {
                        trailingNans = 0;
                    }
                }

                int trimmedSize = xValues.length - leadingNans - trailingNans;
                double[] filteredXValues = getSeriesArray(this.seriesXValues, seriesName, trimmedSize);
                double[] filteredYValues = getSeriesArray(this.seriesYValues, seriesName, trimmedSize);
                for (int j = 0; j < trimmedSize; j++) {
                    double x = xUnit.fromSim(xValues[j + leadingNans]);
                    double y = yUnit.fromSim(data[j + leadingNans]);

                    if (xAxisLog && x <= 0) {
                        x = Double.NaN;
                    }
                    if (yAxisLog && y <= 0) {
                        y = Double.NaN;
                    }
                    filteredXValues[j] = x;
                    filteredYValues[j] = y;
                }

                newXValues[i] = filteredXValues;
                newYValues[i] = filteredYValues;
            }
        }
        SwingUtilities.invokeLater(() -> {
            for (int i = 0; i < nSeries; i++) {
                String name = dataSet.getName(i);
                this.plot.getSeriesMap().get(name).setEnabled(!(newYValues[i].length == 0));
                this.plot.updateXYSeries(name, newXValues[i], newYValues[i], null);
            }
            if (this.layer.isShowing()) {
                this.layer.revalidate();
                this.layer.repaint();
            }
        });
    }

    private static class ZoomUI extends LayerUI<XChartPanel<XYChart>> {
        private final Rectangle2D rect = new Rectangle2D.Double(0, 0, 0, 0);
        private boolean shown = false;

        @Override
        public void paint(Graphics g, JComponent c) {
            super.paint(g, c);

            if (shown) {
                Graphics2D g2 = (Graphics2D) g.create();
                g2.draw(this.rect);
            }
        }

        /**
         * Sets the dimensions of the rectangle to draw
         *
         * @param x x coordinate of the upper left hand corner (x grows to the right)
         * @param y y coordinate of the upper left hand corner (y grows down)
         * @param w the width of the rectangle
         * @param h the height of the rectangle
         */
        public void setRect(int x, int y, int w, int h) {
            rect.setRect(x, y, w, h);
        }

        /**
         * Sets whether the zoom indicator should be shown.
         * @param b true to show, false to hide.
         */
        public void setRectShown(boolean b) {
            shown = b;
        }

        /**
         * Reset the rectangle coordinates to all 0. Should
         * usually be called right after the rectangle is hidden
         * to prevent a flash of the old dimensions once it is shown
         * again.
         */
        public void resetRect() {
            rect.setRect(0, 0, 0, 0);
        }
    }

    private class ZoomMouseAdapter extends MouseAdapter {
        private int pressedX;
        private int pressedY;
        private int currentX;
        private int currentY;
        private boolean drawingRect = false;

        @Override
        public void mouseClicked(MouseEvent e) {
            super.mouseClicked(e);
        }

        @Override
        public void mousePressed(MouseEvent e) {
            if (e.isPopupTrigger()) {
                return;
            }
            this.pressedX = e.getX();
            this.pressedY = e.getY();

            zoomUI.setRect(pressedX, pressedY, 0, 0);
            zoomUI.setRectShown(true);
            this.drawingRect = true;
            layer.repaint();
        }

        @Override
        public void mouseReleased(MouseEvent e) {
            if (e.isPopupTrigger() || !drawingRect) {
                return;
            }

            zoomUI.setRectShown(false);
            this.drawingRect = false;

            if (zoomUI.rect.getWidth() == 0 || zoomUI.rect.getHeight() == 0) {
                layer.repaint();
                zoomUI.resetRect();
                return;
            }
            zoomUI.resetRect();
            int xMin = Math.min(pressedX, currentX);
            int xMax = Math.max(pressedX, currentX);

            int yTop = Math.min(pressedY, currentY);
            int yBottom = Math.max(pressedY, currentY);

            double newXMin = plot.getChartXFromCoordinate(xMin);
            double newXMax = plot.getChartXFromCoordinate(xMax);

            double newYMin = plot.getChartYFromCoordinate(yBottom);
            double newYMax = plot.getChartYFromCoordinate(yTop);

            plot.getStyler().setXAxisMin(newXMin);
            plot.getStyler().setXAxisMax(newXMax);

            plot.getStyler().setYAxisMin(newYMin);
            plot.getStyler().setYAxisMax(newYMax);

            layer.repaint();
        }

        @Override
        public void mouseDragged(MouseEvent e) {
            if (!drawingRect) {
                return;
            }
            currentX = e.getX();
            currentY = e.getY();

            int xDist = Math.abs(currentX - pressedX);
            int yDist = Math.abs(currentY - pressedY);


            zoomUI.setRect(
                    Math.min(pressedX, currentX),
                    Math.min(pressedY, currentY),
                    xDist,
                    yDist
            );
            layer.repaint();
        }
    }
}
