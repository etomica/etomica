/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.heatcapacity;

import etomica.action.IAction;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtom;
import etomica.data.*;
import etomica.data.history.HistoryScrolling;
import etomica.graphics.*;
import etomica.integrator.IntegratorListenerAction;
import etomica.modifier.Modifier;
import etomica.modifier.ModifierGeneral;
import etomica.species.ISpecies;
import etomica.units.Pixel;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Null;
import etomica.units.dimensions.Temperature;
import org.knowm.xchart.CategoryChart;
import org.knowm.xchart.CategoryChartBuilder;
import org.knowm.xchart.XChartPanel;
import org.knowm.xchart.style.Styler;

import javax.swing.*;
import java.awt.*;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.List;
import java.util.function.DoubleConsumer;
import java.util.function.DoubleSupplier;

/**
 * Graphical front end of the Heat Capacity module.
 */
public class HeatCapacityGraphic extends SimulationGraphic {

    private static final String APP_NAME = "Heat Capacity";
    private static final int REPAINT_INTERVAL = 10;
    /** number of integrator steps between points of the temperature plot, one block of HeatCapacityMeter samples */
    private static final int HISTORY_INTERVAL = HeatCapacityMeter.PLOT_BLOCK * HeatCapacityMeter.INTERVAL;
    private static final int HISTORY_LENGTH = 300;
    /** number of integrator steps between checks for updating the metrics and the bar chart */
    private static final int DISPLAY_INTERVAL = 100;
    /** default sleep period (ms), which corresponds to 1 on the delay slider */
    private static final double DEFAULT_SLEEP = 1;
    private static final List<String> MODE_LABELS = Arrays.asList(
            "Translational KE", "Rotational KE", "Vibrational KE", "Vibrational PE", "Intermolecular PE");

    protected final HeatCapacity sim;
    protected final HeatCapacityMeter meter;
    protected final DeviceThermoSlider thermoSlider;
    protected final DeviceSlider nMonoSlider, nDiSlider, nTriSlider;
    protected final JRadioButton monoButton, diButton, triButton;
    protected final ButtonGroup speciesGroup;
    protected final JButton heatButton;
    protected final DeviceSlider heatSlider;
    protected final CategoryChart partitionChart;
    protected final XChartPanel<CategoryChart> partitionPanel;
    protected final JLabel[] metricLabels = new JLabel[7];
    protected final JLabel noteLabel = new JLabel(" ");
    protected final AccumulatorHistory[] histories;
    protected final DataPumpListener[] historyPumps;
    protected double heatAmount = 0.5;
    protected volatile boolean constructed;
    /** molecule counts as of the last time molecules were placed */
    protected int lastMono, lastDi, lastTri;
    /** time that the displays were last refreshed, used to limit the refresh rate */
    protected long lastRefresh;

    public HeatCapacityGraphic(final HeatCapacity simulation) {
        super(simulation, TABBED_PANE, APP_NAME, REPAINT_INTERVAL);
        this.sim = simulation;
        this.meter = sim.meter;

        // slow enough to see energy spread among the modes after adding heat, but the delay slider
        // can remove this for faster data collection.
        sim.getController().setSleepPeriod(DEFAULT_SLEEP);

        // ---- configuration display ----
        DisplayBox displayBox = getDisplayBox(sim.box);
        displayBox.setLabel("Configuration");
        displayBox.setPixelUnit(new Pixel(20));
        // each species has its own colors, and atoms within a molecule differ in color, so that
        // they can be told apart as the molecule rotates and vibrates
        final Color[] monoColors = {new Color(0x2e9e4f)};
        final Color[] diColors = {new Color(0xf2c230), new Color(0xe07020)};
        //final Color[] triColors = {new Color(0xcc3333), Color.WHITE, new Color(0x2a7fff)};
        final Color[] triColors = {Color.WHITE, Color.RED, Color.WHITE};
        displayBox.setColorScheme(new ColorScheme() {
            public Color getAtomColor(IAtom a) {
                ISpecies species = a.getParentGroup().getType();
                Color[] colors = species == sim.speciesTri ? triColors
                        : (species == sim.speciesDi ? diColors : monoColors);
                return colors[a.getIndex() % colors.length];
            }
        });
        DiameterHashByType diameters = new DiameterHashByType();
        diameters.setDiameter(sim.speciesMono.getLeafType(), 1.0);
        diameters.setDiameter(sim.speciesDi.getLeafType(), 1.0);
        diameters.setDiameter(sim.speciesTri.getLeafType(), 1.0);
        displayBox.setDiameterHash(diameters);

        thermoSlider = new DeviceThermoSlider(sim.getController(), sim.integrator);
        final double initialT = sim.integrator.getTemperature();
        thermoSlider.setPrecision(2);
        thermoSlider.setMinimum(0.2);
        thermoSlider.setMaximum(3);
        thermoSlider.setSliderMajorValues(4);
        // changing the scale of the slider sends the integrator a meaningless temperature, perhaps
        // still pending.  Restore the temperature, and a new modifier cancels anything pending.
        sim.integrator.setTemperature(initialT);
        thermoSlider.setModifier(new ModifierGeneral(sim.integrator, "temperature"));
        thermoSlider.setIsothermal();

        // ---- standard buttons ----
        final IAction clearAverages = () -> {
            meter.reset();
            refreshDisplays();
        };
        getController().getSimRestart().setConfiguration(sim.config);
        getController().getSimRestart().setIgnoreOverlap(true);
        getController().getReinitButton().setPostAction(() -> {
            // a lattice start is not an equilibrated adiabatic state, so equilibrate at constant T
            thermoSlider.setIsothermal();
            clearHistories();
            clearAverages.actionPerformed();
            displayBox.repaint();
        });
        getController().getResetAveragesButton().setPostAction(clearAverages);
        final IAction resetData = getController().getSimRestart().getDataResetAction();
        final IAction resetAll = () -> {
            resetData.actionPerformed();
            clearAverages.actionPerformed();
        };

        // ---- State tab ----
        JTabbedPane setupPanel = new JTabbedPane();
        GridBagConstraints vertGBC = SimulationPanel.getVertGBC();

        thermoSlider.setSliderPostAction(resetAll);
        thermoSlider.setRadioGroupPostAction(() -> {
            // the two kinds of measurement apply to data from one mode only
            meter.resetCalorimetry();
            meter.resetFluctuations();
            refreshDisplays();
        });

        heatSlider = makeSlider("Heat per molecule (units of kT)", 0, 2, 2,
                new SimpleModifier("Heat", () -> heatAmount, v -> heatAmount = v));
        heatButton = new JButton("Add heat");
        heatButton.setEnabled(false);
        heatButton.addActionListener(e -> sim.getController().submitActionInterrupt(new IAction() {
            public void actionPerformed() {
                if (!sim.integrator.isIsothermal()) {
                    meter.addHeat(heatAmount);
                    refreshDisplays();
                }
            }
        }));

        JPanel statePanel = new JPanel(new GridBagLayout());
        statePanel.add(thermoSlider.graphic(), vertGBC);
        statePanel.add(heatSlider.graphic(), vertGBC);
        statePanel.add(heatButton, vertGBC);
        setupPanel.add(statePanel, "State");

        // ---- Molecules tab ----
        lastMono = sim.getN(sim.speciesMono);
        lastDi = sim.getN(sim.speciesDi);
        lastTri = sim.getN(sim.speciesTri);
        nMonoSlider = makeCountSlider("Monatomic", sim.speciesMono);
        nDiSlider = makeCountSlider("Diatomic", sim.speciesDi);
        nTriSlider = makeCountSlider("Triatomic", sim.speciesTri);

        monoButton = new JRadioButton("Monatomic");
        diButton = new JRadioButton("Diatomic", true);
        triButton = new JRadioButton("Triatomic");
        speciesGroup = new ButtonGroup();
        speciesGroup.add(monoButton);
        speciesGroup.add(diButton);
        speciesGroup.add(triButton);
        monoButton.addActionListener(e -> selectComposition(HeatCapacity.DEFAULT_MOLECULES, 0, 0));
        diButton.addActionListener(e -> selectComposition(0, HeatCapacity.DEFAULT_MOLECULES, 0));
        triButton.addActionListener(e -> selectComposition(0, 0, HeatCapacity.DEFAULT_MOLECULES));
        JPanel radioPanel = new JPanel(new GridLayout(1, 0));
        radioPanel.add(monoButton);
        radioPanel.add(diButton);
        radioPanel.add(triButton);

        JPanel moleculesPanel = new JPanel(new GridBagLayout());
        moleculesPanel.add(radioPanel, vertGBC);
        moleculesPanel.add(nMonoSlider.graphic(), vertGBC);
        moleculesPanel.add(nDiSlider.graphic(), vertGBC);
        moleculesPanel.add(nTriSlider.graphic(), vertGBC);
        setupPanel.add(moleculesPanel, "Molecules");

        // ---- Potential tab ----
        DeviceSlider kSlider = makeSlider("Bond stiffness", 50, 1000, 0,
                new SimpleModifier("Bond stiffness", sim::getBondStiffness, sim::setBondStiffness));
        kSlider.setPostAction(() -> {
            clearHistories();
            resetAll.actionPerformed();
        });
        JPanel potentialPanel = new JPanel(new GridBagLayout());
        potentialPanel.add(kSlider.graphic(), vertGBC);
        potentialPanel.add(new JLabel("<html><div style='width:220px'>In classical mechanics, the heat capacity "
                + "does not depend on the bond stiffness. Stiffer bonds exchange energy with translation "
                + "more slowly, but they end up with the same share of the energy.</div></html>"), vertGBC);
        setupPanel.add(potentialPanel, "Potential");

        // ---- Controls tab ----
        // DeviceDelaySlider always starts at no delay, so make a similar one that keeps the default
        DeviceSlider delaySlider = new DeviceSlider(sim.getController(), new SimpleModifier("Delay",
                () -> Math.sqrt(sim.getController().getSleepPeriod()), v -> sim.getController().setSleepPeriod(v * v)), 1);
        delaySlider.setShowValues(false);
        delaySlider.setMinimum(0);
        delaySlider.setMaximum(10);
        delaySlider.setNMajor(0);
        delaySlider.setValue(Math.sqrt(DEFAULT_SLEEP));
        delaySlider.setLabel("Simulation Delay");
        delaySlider.setShowBorder(true);
        Hashtable<Integer, JLabel> delayLabels = new Hashtable<>();
        delayLabels.put(0, new JLabel("fast", JLabel.CENTER));
        delayLabels.put(100, new JLabel("slow", JLabel.CENTER)); // slider scales values by 10
        delaySlider.getSlider().setLabelTable(delayLabels);
        JPanel controlsPanel = new JPanel(new GridBagLayout());
        controlsPanel.add(delaySlider.graphic(), vertGBC);
        controlsPanel.add(new JLabel("<html><div style='width:220px'>Move the delay to zero to collect data "
                + "quickly. With a delay, it is easier to follow how energy spreads among the modes after "
                + "adding heat.</div></html>"), vertGBC);
        setupPanel.add(controlsPanel, "Controls");

        getPanel().controlPanel.add(setupPanel, vertGBC);

        // ---- Temperatures tab ----
        DisplayPlotXChart tPlot = new DisplayPlotXChart();
        tPlot.setLabel("Temperatures");
        tPlot.setXLabel("Time");
        tPlot.setYLabel("Temperature");
        DataSourceCountTime timeCounter = new DataSourceCountTime(sim.integrator);
        String[] tLabels = {"T", "T trans", "T rot", "T vib"};
        histories = new AccumulatorHistory[tLabels.length];
        historyPumps = new DataPumpListener[tLabels.length];
        for (int i = 0; i < tLabels.length; i++) {
            final int iT = i;
            HistoryScrolling history = new HistoryScrolling(HISTORY_LENGTH);
            histories[i] = new AccumulatorHistory(history);
            histories[i].setTimeDataSource(timeCounter);
            historyPumps[i] = new DataPumpListener(makeTemperatureSource(tLabels[i], () -> meter.getPlotT(iT)), histories[i], HISTORY_INTERVAL);
            sim.integrator.getEventManager().addListener(historyPumps[i]);
            getController().getDataStreamPumps().add(historyPumps[i]);
            histories[i].setDataSink(tPlot.makeSink(tLabels[i]));
        }
        // XChart hangs if tick labels can't be made distinct, so a pattern must allow for short times
        tPlot.getChart().getStyler().setXAxisDecimalPattern("0.###");
        tPlot.getChart().getStyler().setMarkerSize(0);
        add(tPlot);

        // ---- Energy partition tab ----
        partitionChart = new CategoryChartBuilder().title("Average energy per molecule")
                .xAxisTitle("").yAxisTitle("Energy / kT").build();
        partitionChart.getStyler().setLegendPosition(Styler.LegendPosition.InsideNE);
        partitionChart.getStyler().setYAxisMin(0.0);
        partitionChart.addSeries("Measured", MODE_LABELS, Arrays.asList(0.0, 0.0, 0.0, 0.0, 0.0));
        partitionChart.addSeries("Equipartition", MODE_LABELS, Arrays.asList(0.0, 0.0, 0.0, 0.0, 0.0));
        partitionPanel = new XChartPanel<>(partitionChart);
        addAsTab(partitionPanel, "Energy partition", true);

        // ---- Metrics tab ----
        JPanel metrics = new JPanel(new GridBagLayout());
        GridBagConstraints gbc = SimulationPanel.getVertGBC();
        for (int i = 0; i < metricLabels.length; i++) {
            metricLabels[i] = new JLabel(" ");
            metricLabels[i].setFont(metricLabels[i].getFont().deriveFont(Font.BOLD, 16f));
            metrics.add(metricLabels[i], gbc);
        }
        metrics.add(noteLabel, gbc);
        addAsTab(metrics, "Metrics", true);

        // ---- Reference tab ----
        JEditorPane reference = new JEditorPane("text/html", REFERENCE_TEXT);
        reference.setEditable(false);
        reference.setCaretPosition(0);
        addAsTab(new JScrollPane(reference), "Reference", true);

        // repaint every step when the simulation is slowed, as for the Thermal Conductivity module
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

        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(new IAction() {
            public void actionPerformed() {
                long now = System.currentTimeMillis();
                if (now - lastRefresh < 250) return;
                lastRefresh = now;
                refreshDisplays();
            }
        }, DISPLAY_INTERVAL));
        clearHistories();
        constructed = true;
        refreshDisplays();
    }

    private static final String REFERENCE_TEXT = "<html><body style='font-family:sans-serif;margin:12px'>"
            + "<h3>Equipartition</h3>"
            + "<p>Each quadratic term in the energy contributes kT/2 to the average energy, and so k/2 to the "
            + "heat capacity C<sub>v</sub>.</p>"
            + "<table border='1' cellpadding='4' cellspacing='0'>"
            + "<tr><th>Molecule</th><th>Translation</th><th>Rotation</th><th>Vibration (KE + PE)</th><th>C<sub>v</sub>/k</th></tr>"
            + "<tr><td>Monatomic</td><td align='center'>3/2</td><td align='center'>0</td><td align='center'>0</td><td align='center'><b>3/2</b></td></tr>"
            + "<tr><td>Flexible diatomic</td><td align='center'>3/2</td><td align='center'>1</td><td align='center'>1</td><td align='center'><b>7/2</b></td></tr>"
            + "<tr><td>Flexible bent triatomic</td><td align='center'>3/2</td><td align='center'>3/2</td><td align='center'>3</td><td align='center'><b>6</b></td></tr>"
            + "</table>"
            + "<h3>Experiment, ideal gas at 298 K</h3>"
            + "<table border='1' cellpadding='4' cellspacing='0'>"
            + "<tr><th>Gas</th><th>C<sub>v</sub>/R, experiment</th><th>C<sub>v</sub>/R, classical</th></tr>"
            + "<tr><td>He</td><td align='center'>1.50</td><td align='center'>3/2</td></tr>"
            + "<tr><td>N<sub>2</sub></td><td align='center'>2.50</td><td align='center'>7/2</td></tr>"
            + "<tr><td>H<sub>2</sub>O vapor</td><td align='center'>3.04</td><td align='center'>6</td></tr>"
            + "</table>"
            + "<p>The real gases fall below the classical values because their vibrations are quantum-mechanically "
            + "frozen at ordinary temperatures. A classical simulation such as this one cannot show that.</p>"
            + "<p style='font-size:smaller'>Experimental values from C<sub>p</sub> in the NIST-JANAF tables "
            + "(Chase, 1998), with C<sub>v</sub> = C<sub>p</sub> &minus; R.</p>"
            + "<h3>Two ways to measure C<sub>v</sub></h3>"
            + "<p><b>Calorimetric</b> (adiabatic): C<sub>v</sub> = Q / (N &Delta;T), using the rise in the average "
            + "temperature after adding heat Q.</p>"
            + "<p><b>Fluctuation</b> (isothermal): C<sub>v</sub> = (&lt;E<sup>2</sup>&gt; &minus; &lt;E&gt;<sup>2</sup>) / "
            + "(N k T<sup>2</sup>), using the fluctuations of the total energy E.</p>"
            + "<h3>Using the calorimeter</h3>"
            + "<p>Select Adiabatic, let the simulation run for a little while to get the starting temperature, then "
            + "press Add heat. The heat goes only into translation. Watch the Temperatures tab as "
            + "collisions share it with rotation and vibration.</p>"
            + "</body></html>";

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

    protected DeviceSlider makeCountSlider(String label, final ISpecies species) {
        DeviceSlider slider = makeSlider(label, 0, 300, 0,
                new SimpleModifier(label, () -> sim.getN(species), v -> sim.setN(species, (int) Math.round(v))));
        slider.setPostAction(this::compositionChanged);
        return slider;
    }

    /**
     * Selects a composition using the radio buttons.  Setting the sliders changes the simulation
     * (via their modifiers), one species at a time.
     */
    protected void selectComposition(int nMono, int nDi, int nTri) {
        nMonoSlider.setValue(nMono);
        nDiSlider.setValue(nDi);
        nTriSlider.setValue(nTri);
    }

    /**
     * Responds to changes from the molecule sliders, placing molecules anew (which also resets the
     * averages and histories) and returning to isothermal, because a lattice start is not an
     * equilibrated adiabatic state.  Also updates the radio buttons, none of which is
     * selected for a mixture.
     */
    protected void compositionChanged() {
        final int nm = sim.getN(sim.speciesMono), nd = sim.getN(sim.speciesDi), nt = sim.getN(sim.speciesTri);
        if (nm == lastMono && nd == lastDi && nt == lastTri) return;
        lastMono = nm;
        lastDi = nd;
        lastTri = nt;
        getController().getSimRestart().actionPerformed();
        thermoSlider.setIsothermal();
        clearHistories();
        meter.reset();
        refreshDisplays();
        getDisplayBox(sim.box).repaint();
        final int nDefault = HeatCapacity.DEFAULT_MOLECULES;
        SwingUtilities.invokeLater(() -> {
            if (nm == nDefault && nd == 0 && nt == 0) monoButton.setSelected(true);
            else if (nm == 0 && nd == nDefault && nt == 0) diButton.setSelected(true);
            else if (nm == 0 && nd == 0 && nt == nDefault) triButton.setSelected(true);
            else speciesGroup.clearSelection();
        });
    }

    /**
     * Starts the temperature plot over, leaving a point for each temperature because the plot
     * can't paint if all of its series are empty.
     */
    protected void clearHistories() {
        meter.restartPlot();
        for (int i = 0; i < histories.length; i++) {
            histories[i].reset();
            historyPumps[i].actionPerformed();
        }
    }

    protected DataSourceScalar makeTemperatureSource(String label, final DoubleSupplier getter) {
        return new DataSourceScalar(label, Temperature.DIMENSION) {
            public double getDataAsScalar() {
                return getter.getAsDouble();
            }
        };
    }

    /**
     * Gathers the current measurements (called from the integrator/controller thread) and
     * schedules the display update on the Swing thread.
     */
    protected void refreshDisplays() {
        if (!constructed) return; // actions pending from slider setup might get here first
        final boolean isothermal = sim.integrator.isIsothermal();
        final double T = meter.getTemperature();
        final double avgT = meter.getAverageT();
        final double[] cal = meter.getCalorimetricCv();
        final double[] fluct = meter.getFluctuationCv();
        final double pred = meter.getPredictedCv();
        final double[] partition = meter.getPartition();
        final double[] equipartition = meter.getEquipartition();
        final boolean heatAdded = meter.isHeatAdded(), hasBaseline = meter.hasBaseline();
        final double sinceHeat = meter.getTimeSinceHeat();
        final int nm = sim.getN(sim.speciesMono), nd = sim.getN(sim.speciesDi), nt = sim.getN(sim.speciesTri);
        SwingUtilities.invokeLater(() -> {
            heatButton.setEnabled(!isothermal && nm + nd + nt > 1);
            metricLabels[0].setText("T (current) = " + fmt(T, 3));
            metricLabels[1].setText("T (average since reset or heat addition) = " + fmt(avgT, 3));
            metricLabels[2].setText("Cv/k per molecule, calorimetric = " + fmtErr(cal));
            metricLabels[3].setText("Cv/k per molecule, fluctuation = " + (isothermal ? fmtErr(fluct) : "isothermal only"));
            metricLabels[4].setText("Cv/k per molecule, equipartition = " + fmt(pred, 3));
            metricLabels[5].setText("Monatomic: " + nm + "    Diatomic: " + nd + "    Triatomic: " + nt);
            metricLabels[6].setText(" ");
            noteLabel.setText("<html>" + calorimetryNote(isothermal, heatAdded, hasBaseline, sinceHeat) + "</html>");
            partitionChart.updateCategorySeries("Measured", MODE_LABELS, finiteList(partition), null);
            partitionChart.updateCategorySeries("Equipartition", MODE_LABELS, finiteList(equipartition), null);
            partitionPanel.repaint();
        });
    }

    private static String calorimetryNote(boolean isothermal, boolean heatAdded, boolean hasBaseline, double sinceHeat) {
        if (isothermal) {
            return "To measure Cv with a calorimeter, select Adiabatic and then add heat.";
        }
        if (!heatAdded) {
            return "Adiabatic: the starting temperature is being measured. Let this run a little, then press Add heat.";
        }
        if (!hasBaseline) {
            return "Heat was added before the starting temperature was measured well enough to give Cv.<br>"
                    + "Let this run a little, and add heat again.";
        }
        if (0.5 * sinceHeat < HeatCapacityMeter.MIN_WINDOW) {
            return "Heat added. Waiting for the energy to spread among the modes.";
        }
        return "Calorimetric Cv uses the average T over the second half of the time since heat was added,<br>"
                + "so the estimate improves as the run continues.";
    }

    private static List<Double> finiteList(double[] v) {
        Double[] d = new Double[v.length];
        for (int i = 0; i < v.length; i++) d[i] = Double.isNaN(v[i]) ? 0.0 : v[i];
        return Arrays.asList(d);
    }

    protected static String fmt(double v, int digits) {
        return Double.isNaN(v) || Double.isInfinite(v) ? "—" : String.format("%." + digits + "f", v);
    }

    protected static String fmtErr(double[] v) {
        if (Double.isNaN(v[0])) return "—";
        return fmt(v[0], 2) + (Double.isNaN(v[1]) ? "" : " ± " + fmt(v[1], 2));
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
        HeatCapacityGraphic graphic = new HeatCapacityGraphic(new HeatCapacity());
        SimulationGraphic.makeAndDisplayFrame(graphic.getPanel(), APP_NAME);
    }
}
