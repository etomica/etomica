package etomica.graphics;

import org.knowm.xchart.XYChart;
import org.knowm.xchart.XYSeries;
import org.knowm.xchart.internal.chartpart.ToolTips;

import javax.swing.*;
import javax.swing.event.PopupMenuEvent;
import javax.swing.event.PopupMenuListener;
import java.awt.*;
import java.awt.event.*;
import java.util.Map;

class ChartPanel extends JPanel {
    private final XYChart chart;
    private final DisplayPlotXChart display;
    private int oldMarkerSize;

    ChartPanel(DisplayPlotXChart display) {
        this.display = display;
        this.chart = display.getChart();
        this.oldMarkerSize = this.chart.getStyler().getMarkerSize();
        this.setPreferredSize(new Dimension(chart.getWidth(), chart.getHeight()));

        ToolTips toolTips = chart.getToolTips();
        if (toolTips != null) {
            MouseMotionListener mml = toolTips.getMouseMotionListener();
            if (mml != null) {
                this.addMouseMotionListener(mml);
            }
        }

        this.addMouseMotionListener(chart.getCursor());

        this.addMouseWheelListener(e -> {
            if (!e.isShiftDown()) {
                e.getComponent().getParent().dispatchEvent(e);
                return;
            }
            double scale = e.getPreciseWheelRotation() * -0.05;
            System.out.println(scale);

            double chartX = e.getX();
            double chartY = e.getY();
            double chartXMin = chart.getStyler().getXAxisMin() == null ? display.getDataXMin() : chart.getStyler().getXAxisMin();
            double chartYMin = chart.getStyler().getYAxisMin() == null ? display.getDataYMin() : chart.getStyler().getYAxisMin();
            double chartXMax = chart.getStyler().getXAxisMax() == null ? display.getDataXMax() : chart.getStyler().getXAxisMax();
            double chartYMax = chart.getStyler().getYAxisMax() == null ? display.getDataYMax() : chart.getStyler().getYAxisMax();

            if (Double.isNaN(chartXMin)) {
                return;
            }

            chartXMin = chart.getScreenXFromChart(chartXMin);
            chartYMin = chart.getScreenYFromChart(chartYMin);
            chartXMax = chart.getScreenXFromChart(chartXMax);
            chartYMax = chart.getScreenYFromChart(chartYMax);

            double xRange = (chartXMax - chartXMin);
            double yRange = (chartYMax - chartYMin);

            double xFrac = (chartX - chartXMin) / xRange;
            double yFrac = (chartY - chartYMin) / yRange;

            chartXMin += scale * xRange * xFrac;
            chartYMin += scale * yRange * yFrac;
            chartXMax -= scale * xRange * (1 - xFrac);
            chartYMax -= scale * yRange * (1 - yFrac);

            chartXMin = chart.getChartXFromCoordinate((int) chartXMin);
            chartYMin = chart.getChartYFromCoordinate((int) chartYMin);
            chartXMax = chart.getChartXFromCoordinate((int) chartXMax);
            chartYMax = chart.getChartYFromCoordinate((int) chartYMax);

            chart.getStyler().setXAxisMin(chartXMin);
            chart.getStyler().setXAxisMax(chartXMax);
            chart.getStyler().setYAxisMin(chartYMin);
            chart.getStyler().setYAxisMax(chartYMax);
            display.getPanel().repaint();
        });

        this.addMouseListener(new MouseAdapter() {
            @Override
            public void mousePressed(MouseEvent e) {
                if (e.isPopupTrigger()) {
                    showMenu(e);
                }
            }

            @Override
            public void mouseReleased(MouseEvent e) {
                if (e.isPopupTrigger()) {
                    showMenu(e);
                }
            }

            private void showMenu(MouseEvent e) {
                ChartContextMenu menu = new ChartContextMenu();
                menu.show(e.getComponent(), e.getX(), e.getY());
                menu.getGraphics().dispose();
            }
        });

    }


    @Override
    protected void paintComponent(Graphics g) {
        super.paintComponent(g);

        Graphics2D g2d = (Graphics2D) g.create();
        this.chart.paint(g2d, this.getWidth(), this.getHeight());
        g2d.dispose();
    }

    private class ChartContextMenu extends JPopupMenu {
        public ChartContextMenu() {
            JMenuItem resetZoomItem = new JMenuItem("Reset Zoom");
            resetZoomItem.addActionListener(e -> {
                display.resetZoom();
            });
            this.add(resetZoomItem);

            JCheckBoxMenuItem xLogItem = new JCheckBoxMenuItem("X Axis Log", chart.getStyler().isXAxisLogarithmic());
            xLogItem.addActionListener(e -> {
                display.resetZoom();
                display.clearData();
                chart.getStyler().setXAxisLogarithmic(!chart.getStyler().isXAxisLogarithmic());
                display.doUpdate();
                xLogItem.setState(chart.getStyler().isXAxisLogarithmic());
            });
            this.add(xLogItem);

            JCheckBoxMenuItem yLogItem = new JCheckBoxMenuItem("Y Axis Log", chart.getStyler().isYAxisLogarithmic());
            yLogItem.addActionListener(e -> {
                display.resetZoom();
                display.clearData();
                chart.getStyler().setYAxisLogarithmic(!chart.getStyler().isYAxisLogarithmic());
                display.doUpdate();
                yLogItem.setState(chart.getStyler().isYAxisLogarithmic());
            });
            this.add(yLogItem);

            JCheckBoxMenuItem togglePointsItem = new JCheckBoxMenuItem("Display Points", chart.getStyler().getMarkerSize() != 0);
            togglePointsItem.addActionListener(e -> {
                if (chart.getStyler().getMarkerSize() == 0) {
                    chart.getStyler().setMarkerSize(oldMarkerSize);
                } else {
                    oldMarkerSize = chart.getStyler().getMarkerSize();
                    chart.getStyler().setMarkerSize(0);
                }
                togglePointsItem.setState(chart.getStyler().getMarkerSize() != 0);
                ChartPanel.this.repaint();
            });
            this.add(togglePointsItem);

            this.add(new JSeparator());

            JMenu displayDataMenu = new JMenu("Show data");
            // TODO select all
            this.addPopupMenuListener(new PopupMenuListener() {
                @Override
                public void popupMenuWillBecomeVisible(PopupMenuEvent popupMenuEvent) {
                    displayDataMenu.removeAll();
                    for (Map.Entry<String, XYSeries> entry : chart.getSeriesMap().entrySet()) {
                        XYSeries series = entry.getValue();
                        String label = series.getLabel();
                        JMenuItem menuItem = new JMenuItem(label);
                        menuItem.addActionListener(e -> {
                            JFrame frame = new JFrame();
                            JTextArea textArea = new JTextArea();
                            JScrollPane scrollPane = new JScrollPane(textArea);
                            textArea.setEditable(false);

                            textArea.append(chart.getXAxisTitle() + "\t" + label + "\n");

                            double[] xValues = series.getXData();
                            double[] yValues = series.getYData();

                            for (int i = 0; i < xValues.length; i++) {
                                textArea.append(xValues[i] + "\t" + yValues[i] + "\n");
                            }
                            frame.add(scrollPane);
                            frame.pack();
                            frame.setVisible(true);
                        });
                        displayDataMenu.add(menuItem);
                    }
                }

                @Override
                public void popupMenuWillBecomeInvisible(PopupMenuEvent e) {

                }

                @Override
                public void popupMenuCanceled(PopupMenuEvent e) {

                }
            });

            this.add(displayDataMenu);


        }
    }
}
