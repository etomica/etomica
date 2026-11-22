/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.thermalconductivity;

import etomica.action.IAction;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.graphics.*;
import etomica.integrator.IntegratorListenerAction;
import etomica.modifier.Modifier;
import etomica.units.Pixel;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Null;
import org.knowm.xchart.XChartPanel;
import org.knowm.xchart.XYChart;
import org.knowm.xchart.XYChartBuilder;
import org.knowm.xchart.XYSeries;
import org.knowm.xchart.style.Styler;
import org.knowm.xchart.style.markers.SeriesMarkers;

import javax.swing.*;
import java.awt.*;
import java.util.ArrayList;
import java.util.List;
import java.util.function.DoubleConsumer;
import java.util.function.DoubleSupplier;

/**
 * Graphical front end of the Thermal Conductivity module.
 */
public class ThermalConductivityGraphic extends SimulationGraphic {

    private static final String APP_NAME = "Thermal Conductivity";
    private static final int REPAINT_INTERVAL = 10;
    private static final int DATA_INTERVAL = 20;
    /** dilute-gas Chapman-Enskog value of k/sqrt(T) for hard disks (sigma = m = kB = 1): 2/sqrt(pi) */
    private static final double K0_OVER_SQRT_T = 2.0 / Math.sqrt(Math.PI);

    protected final ThermalConductivity sim;
    protected final ConductionMeter meter;

    protected XYChart profileChart, fluxChart, kChart;
    protected final JLabel[] metricLabels = new JLabel[8];
    protected final JLabel noteLabel = new JLabel(" ");
    protected DeviceButton recordButton;
    protected JComboBox<String> fluidChooser;
    protected final List<double[]> recorded = new ArrayList<>(); // {phi, k/sqrtT, err/sqrtT}

    public ThermalConductivityGraphic(final ThermalConductivity simulation) {
        super(simulation, TABBED_PANE, APP_NAME, REPAINT_INTERVAL);
        this.sim = simulation;
        this.meter = sim.meter;

        sim.getController().setSleepPeriod(1);

        // configuration display: disks colored by kinetic energy, hot and cold walls marked
        DisplayBox displayBox = getDisplayBox(sim.box);
        displayBox.setLabel("Configuration");
        displayBox.setPixelUnit(new Pixel(15));
        ((DiameterHashByType) displayBox.getDiameterHash()).setDiameter(sim.species.getLeafType(), ThermalConductivity.SIGMA);
        displayBox.setColorScheme(new KineticEnergyColorScheme());
        displayBox.addDrawable(new WallDrawable());

        // Reinitialize button uses our configuration; both buttons refresh the displays
        getTrioControllerButton().getSimRestart().setConfiguration(sim.makeConfiguration());
        final IAction refresh = this::refreshDisplays;
        getController().getReinitButton().setPostAction(refresh);
        getController().getResetAveragesButton().setPostAction(new IAction() {
            public void actionPerformed() {
                meter.reset();
                refreshDisplays();
            }
        });

        // ---- control panel ----
        JTabbedPane setupPanel = new JTabbedPane();

        DeviceSlider tHotSlider = makeSlider("T hot", 0.1, 5, 2, new SimpleModifier("T hot", sim.wall::getTHot,
                v -> {
                    sim.setWallTemperatures(v, sim.wall.getTCold());
                    refreshDisplays();
                }));
        DeviceSlider tColdSlider = makeSlider("T cold", 0.1, 5, 2, new SimpleModifier("T cold", sim.wall::getTCold,
                v -> {
                    sim.setWallTemperatures(sim.wall.getTHot(), v);
                    refreshDisplays();
                }));
        DeviceSlider phiSlider = makeSlider("Area fraction φ", 0.02, 0.60, 2, new SimpleModifier("phi", sim::getPhi, v -> {
            sim.setPhi(v);
            refreshDisplays();
        }));
        phiSlider.setNMajor(4);

        JPanel statePanel = new JPanel(new GridBagLayout());
        GridBagConstraints vertGBC = SimulationPanel.getVertGBC();
        statePanel.add(tHotSlider.graphic(), vertGBC);
        statePanel.add(tColdSlider.graphic(), vertGBC);
        statePanel.add(phiSlider.graphic(), vertGBC);
        setupPanel.add(statePanel, "State");

        fluidChooser = new JComboBox<>(new String[]{"Ideal gas", "Hard disk"});
        fluidChooser.setSelectedItem("Hard disk");
        fluidChooser.addActionListener(e -> {
            final boolean ideal = "Ideal gas".equals(fluidChooser.getSelectedItem());
            sim.getController().submitActionInterrupt(new IAction() {
                public void actionPerformed() {
                    sim.setIdealGas(ideal);
                    refreshDisplays();
                }
            });
        });
        JPanel potentialPanel = new JPanel(new GridBagLayout());
        potentialPanel.add(fluidChooser, vertGBC);
        setupPanel.add(potentialPanel, "Potential");

        DeviceDelaySlider delaySlider = new DeviceDelaySlider(sim.getController());
        delaySlider.setDelayExponent(3);
        delaySlider.setMaxSleep(100);
        recordButton = new DeviceButton(sim.getController(), new IAction() {
            public void actionPerformed() {
                recordK();
            }
        });
        recordButton.setLabel("Record k");
        DeviceButton clearButton = new DeviceButton(sim.getController(), new IAction() {
            public void actionPerformed() {
                recorded.clear();
                refreshKPlot();
            }
        });
        clearButton.setLabel("Clear recorded points");
        JPanel controlsPanel = new JPanel(new GridBagLayout());
        controlsPanel.add(delaySlider.graphic(), vertGBC);
        controlsPanel.add(recordButton.graphic(), vertGBC);
        controlsPanel.add(clearButton.graphic(), vertGBC);
        setupPanel.add(controlsPanel, "Controls");

        getPanel().controlPanel.add(setupPanel, vertGBC);

        // ---- graphics tabs ----
        profileChart = makeChart("Temperature profile", "x", "T");
        addLine(profileChart, "T(x)", Color.BLACK, true, false);
        addLine(profileChart, "linear fit", new Color(0x2a7fff), false, true);
        addLine(profileChart, "walls", new Color(0xcc3333), true, false);
        profileChart.getSeriesMap().get("walls").setMarkerColor(new Color(0xcc3333));
        profileChart.getStyler().setXAxisMin(0.0).setXAxisMax(sim.getLx());
        addAsTab(new XChartPanel<>(profileChart), "Profile", true);

        fluxChart = makeChart("Heat flux", "time", "J");
        addLine(fluxChart, "J hot", new Color(0xcc3333), false, true);
        addLine(fluxChart, "J cold", new Color(0x2a7fff), false, true);
        addAsTab(new XChartPanel<>(fluxChart), "Flux", true);

        kChart = makeChart("k vs density", "area fraction φ", "k / √T");
        XYSeries kSeries = kChart.addSeries("recorded", new double[]{0}, new double[]{0}, new double[]{0}); // finite placeholders: NaN error bars break XChart's axis ticks
        kSeries.setXYSeriesRenderStyle(XYSeries.XYSeriesRenderStyle.Scatter);
        kSeries.setMarker(SeriesMarkers.CIRCLE);
        kSeries.setMarkerColor(Color.BLACK);
        kSeries.setEnabled(false);
        XYSeries ce = kChart.addSeries("dilute-gas Chapman-Enskog", new double[]{0, 0.65}, new double[]{K0_OVER_SQRT_T, K0_OVER_SQRT_T});
        ce.setXYSeriesRenderStyle(XYSeries.XYSeriesRenderStyle.Line);
        ce.setMarker(SeriesMarkers.NONE);
        ce.setLineColor(new Color(0x888888));
        kChart.getStyler().setXAxisMin(0.0).setXAxisMax(0.65).setYAxisMin(0.0);
        addAsTab(new XChartPanel<>(kChart), "k vs density", true);

        JPanel metrics = new JPanel(new GridBagLayout());
        GridBagConstraints gbc = SimulationPanel.getVertGBC();
        for (int i = 0; i < metricLabels.length; i++) {
            metricLabels[i] = new JLabel(" ");
            metricLabels[i].setFont(metricLabels[i].getFont().deriveFont(Font.BOLD, 16f));
            metrics.add(metricLabels[i], gbc);
        }
        metrics.add(noteLabel, gbc);
        addAsTab(metrics, "Metrics", true);

        // Smooth animation: when the user slows the simulation (delay > 0) repaint every step so
        // disks move a small fraction of a diameter between frames; with no delay repaint only
        // every REPAINT_INTERVAL steps so that equilibration stays fast.
        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(new IAction() {
            private int current = REPAINT_INTERVAL;

            public void actionPerformed() {
                int want = sim.getController().getSleepPeriod() > 0 ? 1 : REPAINT_INTERVAL;
                if (want != current) {
                    current = want;
                    setPaintInterval(sim.box, want);
                }
            }
        }));

        // periodic data update (on the integrator thread; Swing updates are deferred to the EDT)
        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(new IAction() {
            public void actionPerformed() {
                refreshDisplays();
            }
        }, DATA_INTERVAL));
        refreshDisplays();
    }

    protected DeviceSlider makeSlider(String label, double min, double max, int precision, Modifier modifier) {
        DeviceSlider slider = new DeviceSlider(sim.getController(), modifier, precision);
        slider.setMinimum(min);
        slider.setMaximum(max);
        slider.setLabel(label);
        slider.setShowBorder(true);
        slider.setShowValues(true);
        slider.setEditValues(true);
        slider.setNMajor(4);
        return slider;
    }

    protected XYChart makeChart(String title, String xLabel, String yLabel) {
        XYChart chart = new XYChartBuilder().title(title).xAxisTitle(xLabel).yAxisTitle(yLabel).build();
        chart.getStyler().setLegendPosition(Styler.LegendPosition.InsideNE);
        return chart;
    }

    protected void addLine(XYChart chart, String name, Color color, boolean markers, boolean line) {
        XYSeries s = chart.addSeries(name, new double[]{Double.NaN}, new double[]{Double.NaN});
        s.setXYSeriesRenderStyle(line ? XYSeries.XYSeriesRenderStyle.Line : XYSeries.XYSeriesRenderStyle.Scatter);
        s.setMarker(markers ? SeriesMarkers.CIRCLE : SeriesMarkers.NONE);
        s.setLineColor(color);
        s.setMarkerColor(color);
        s.setEnabled(false);
    }

    /** Updates a series with the finite points of (x,y); disables it if there are none. EDT only. */
    protected static void setSeries(XYChart chart, String name, double[] x, double[] y, double[] err) {
        int n = 0;
        for (int i = 0; i < x.length; i++) if (!Double.isNaN(x[i]) && !Double.isNaN(y[i])) n++;
        XYSeries s = chart.getSeriesMap().get(name);
        if (n == 0) {
            s.setEnabled(false);
            return;
        }
        double[] xs = new double[n], ys = new double[n], es = err == null ? null : new double[n];
        int j = 0;
        for (int i = 0; i < x.length; i++) {
            if (Double.isNaN(x[i]) || Double.isNaN(y[i])) continue;
            xs[j] = x[i];
            ys[j] = y[i];
            if (es != null) es[j] = Double.isNaN(err[i]) ? 0 : err[i];
            j++;
        }
        s.setEnabled(true);
        chart.updateXYSeries(name, xs, ys, es);
    }

    /**
     * Gathers the current measurements (called from the integrator/controller thread) and
     * schedules the display update on the Swing thread.
     */
    protected void refreshDisplays() {
        final double[] xb = meter.binCenters();
        final double[] profile = meter.getProfile();
        final double[] fit = meter.getFit();
        final double[] kk = meter.getConductivity();
        final double[][] blocks = meter.getBlocks();
        final double jHot = meter.getJHot(), jCold = meter.getJCold();
        final double time = meter.getElapsedTime();
        final double tHot = sim.wall.getTHot(), tCold = sim.wall.getTCold();
        final double lx = sim.getLx();
        final boolean ideal = sim.isIdealGas();
        final double phi = sim.getPhi();
        final int n = sim.box.getLeafList().size();
        SwingUtilities.invokeLater(() -> {
            setSeries(profileChart, "T(x)", xb, profile, null);
            if (fit != null && !ideal) {
                setSeries(profileChart, "linear fit", new double[]{0, lx},
                        new double[]{fit[1], fit[1] + fit[0] * lx}, null);
            } else {
                setSeries(profileChart, "linear fit", new double[0], new double[0], null);
            }
            setSeries(profileChart, "walls", new double[]{0, lx}, new double[]{tHot, tCold}, null);
            profileChart.getStyler().setYAxisMin(0.0).setYAxisMax(1.1 * Math.max(tHot, tCold));

            double[] bt = new double[blocks.length], bh = new double[blocks.length], bc = new double[blocks.length];
            for (int i = 0; i < blocks.length; i++) {
                bt[i] = blocks[i][0];
                bh[i] = blocks[i][1];
                bc[i] = blocks[i][2];
            }
            setSeries(fluxChart, "J hot", bt, bh, null);
            setSeries(fluxChart, "J cold", bt, bc, null);

            String dash = "—";
            metricLabels[0].setText("J hot = " + fmt(jHot, 4));
            metricLabels[1].setText("J cold = " + fmt(jCold, 4));
            if (ideal) {
                metricLabels[2].setText("dT/dx: no Fourier gradient");
                metricLabels[3].setText("Conductance J/(T hot − T cold) = " + fmt(0.5 * (jHot + jCold) / (tHot - tCold), 4));
                noteLabel.setText("<html>Ideal gas: no collisions, so no interior temperature gradient develops.<br>"
                        + "Molecules carry energy ballistically wall to wall; the flux is set by the wall temperatures.</html>");
            } else {
                metricLabels[2].setText("dT/dx = " + (fit == null ? dash : fmt(fit[0], 5)));
                String kText = Double.isNaN(kk[0]) ? dash : fmt(kk[0], 3) + (Double.isNaN(kk[1]) ? "" : " ± " + fmt(kk[1], 3));
                metricLabels[3].setText("k = " + kText);
                noteLabel.setText("<html>k is computed from the interior temperature gradient (excluding the bins at each wall),<br>"
                        + "not from (T hot − T cold)/L, because of the temperature jump at the walls.</html>");
            }
            metricLabels[4].setText("φ = " + fmt(phi, 3));
            metricLabels[5].setText("N = " + n);
            metricLabels[6].setText("elapsed time = " + fmt(time, 1));
            metricLabels[7].setText(" ");
            recordButton.getButton().setEnabled(!ideal);
            profileChart.setTitle(ideal ? "Temperature profile (ideal gas)" : "Temperature profile");
            // repaint the charts' panels
            for (XYChart c : new XYChart[]{profileChart, fluxChart, kChart}) {
                Container parent = findPanel(c);
                if (parent != null) parent.repaint();
            }
        });
    }

    private final java.util.Map<XYChart, Container> chartPanels = new java.util.HashMap<>();

    private Container findPanel(XYChart chart) {
        Container c = chartPanels.get(chart);
        if (c == null) {
            Component[] comps = getPanel().tabbedPane.getComponents();
            for (Component comp : comps) {
                if (comp instanceof XChartPanel && ((XChartPanel<?>) comp).getChart() == chart) {
                    c = (Container) comp;
                    chartPanels.put(chart, c);
                    break;
                }
            }
        }
        return c;
    }

    protected static String fmt(double v, int digits) {
        return Double.isNaN(v) || Double.isInfinite(v) ? "—" : String.format("%." + digits + "f", v);
    }

    protected void recordK() {
        if (sim.isIdealGas()) return;
        double[] kk = meter.getConductivity();
        if (Double.isNaN(kk[0])) return;
        double sqrtT = Math.sqrt(kk[2]);
        recorded.add(new double[]{sim.getPhi(), kk[0] / sqrtT, kk[1] / sqrtT});
        refreshKPlot();
    }

    protected void refreshKPlot() {
        final double[][] pts = recorded.toArray(new double[0][]);
        SwingUtilities.invokeLater(() -> {
            double[] x = new double[pts.length], y = new double[pts.length], e = new double[pts.length];
            for (int i = 0; i < pts.length; i++) {
                x[i] = pts[i][0];
                y[i] = pts[i][1];
                e[i] = pts[i][2];
            }
            setSeries(kChart, "recorded", x, y, e);
            Container p = findPanel(kChart);
            if (p != null) p.repaint();
        });
    }

    /** Colors disks from blue (low kinetic energy) to red (high). */
    protected class KineticEnergyColorScheme extends ColorScheme {
        public Color getAtomColor(IAtom a) {
            double ke = 0.5 * a.getType().getMass() * ((IAtomKinetic) a).getVelocity().squared();
            double keMax = 2.0 * Math.max(sim.wall.getTHot(), sim.wall.getTCold());
            float f = (float) Math.max(0, Math.min(1, ke / keMax));
            return new Color(f, 0f, 1f - f);
        }
    }

    /** Marks the hot (left, red) and cold (right, blue) walls on the box display. */
    protected class WallDrawable implements Drawable {
        public void draw(Graphics g, int[] origin, double toPixels) {
            int h = (int) (sim.getLy() * toPixels);
            int w = (int) (sim.getLx() * toPixels);
            g.setColor(new Color(0xcc3333));
            g.fillRect(origin[0] - 3, origin[1], 4, h);
            g.setColor(new Color(0x2a7fff));
            g.fillRect(origin[0] + w - 1, origin[1], 4, h);
        }
    }

    /** Modifier built from a getter and setter, for reduced (dimensionless) quantities. */
    protected static class SimpleModifier implements Modifier {
        private final String label;
        private final DoubleSupplier getter;
        private final DoubleConsumer setter;

        SimpleModifier(String label, DoubleSupplier getter, DoubleConsumer setter) {
            this.label = label;
            this.getter = getter;
            this.setter = setter;
        }

        public void setValue(double newValue) {
            setter.accept(newValue);
        }

        public double getValue() {
            return getter.getAsDouble();
        }

        public Dimension getDimension() {
            return Null.DIMENSION;
        }

        public String getLabel() {
            return label;
        }
    }

    public static void main(String[] args) {
        ThermalConductivityGraphic graphic = new ThermalConductivityGraphic(new ThermalConductivity());
        SimulationGraphic.makeAndDisplayFrame(graphic.getPanel(), APP_NAME);
    }
}
