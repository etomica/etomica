/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.nucleation;

import etomica.action.BoxImposePbc;
import etomica.action.IAction;
import etomica.atom.AtomTest;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtom;
import etomica.config.ConfigurationLattice;
import etomica.data.*;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterPressureHard;
import etomica.data.meter.MeterTemperature;
import etomica.data.types.DataDouble;
import etomica.exception.ConfigurationOverlapException;
import etomica.graphics.*;
import etomica.integrator.IntegratorMD;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.modifier.Modifier;
import etomica.modifier.ModifierBoolean;
import etomica.modules.swmd.Swmd;
import etomica.potential.P2SquareWell;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space2d.Space2D;
import etomica.space3d.Space3D;
import etomica.units.*;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Length;
import etomica.units.dimensions.Null;
import etomica.units.dimensions.Pressure2D;
import etomica.util.Constants.CompassDirection;
import net.miginfocom.swing.MigLayout;

import javax.swing.*;
import java.awt.*;
import java.util.ArrayList;

public class NucleationGraphic extends SimulationGraphic {

    private final static String APP_NAME = "Nucleation";
    private final static int REPAINT_INTERVAL = 1;
    protected DeviceThermoSlider tempSlider;
    protected Unit tUnit, eUnit, dUnit, pUnit;
    protected Swmd sim;
    protected final MeterClusterSizes meterClusterSizes;
    protected final MeterLargestCluster meterLargestCluster;

    public NucleationGraphic(final Swmd simulation) {

        super(simulation, TABBED_PANE, APP_NAME, REPAINT_INTERVAL);

        ArrayList<DataPump> dataStreamPumps = getController().getDataStreamPumps();

        final IAction resetDataAction = new IAction() {
            public void actionPerformed() {
                getController().getResetAveragesButton().press();
            }
        };

        this.sim = simulation;

        tUnit = Kelvin.UNIT;

        eUnit = tUnit;
        if (sim.getSpace().D() == 2) {
            dUnit = new SimpleUnit(Null.DIMENSION, 0.344, "Reduced Density", "x", false);
            pUnit = new SimpleUnit(Pressure2D.DIMENSION, 0.173, "Reduced Pressure", "P/Pc", false);
        } else {
            dUnit = new UnitRatio(Mole.UNIT, Liter.UNIT);
            pUnit = Bar.UNIT;
        }

        ((P2SquareWell) sim.potentialWrapper.getWrappedPotential()).setEpsilon(tUnit.toSim(600));
        ((P2SquareWell) sim.potentialWrapper.getWrappedPotential()).setLambda(1.5);

        if (sim.getSpace().D() == 2) {
            int N = 400;
            sim.box.setNMolecules(sim.species, N);
            sim.box.getBoundary().setBoxSize(Vector.of(100, 100));
            new ConfigurationLattice(new LatticeOrthorhombicHexagonal(space), space).initializeCoordinates(sim.box);
            getDisplayBox(sim.box).setPixelUnit(new Pixel(400 / sim.box.getBoundary().getBoxSize().getX(1)));
        } else {
            getDisplayBox(sim.box).setPixelUnit(new Pixel(40 / sim.box.getBoundary().getBoxSize().getX(1)));
        }

        ColorScheme colorSchemeRed = getDisplayBox(sim.box).getColorScheme();
        ColorSchemeCluster colorSchemeCluster = new ColorSchemeCluster(sim.box, new AtomTest() {
            @Override
            public boolean test(IAtom a) {
                return true;
            }
        });
        DeviceCheckBox colorCheckbox = new DeviceCheckBox("color clusters", new ModifierBoolean() {
            @Override
            public void setBoolean(boolean b) {
                getDisplayBox(sim.box).setColorScheme(b ? colorSchemeCluster : colorSchemeRed);
            }

            @Override
            public boolean getBoolean() {
                return getDisplayBox(sim.box).getColorScheme() == colorSchemeCluster;
            }
        });
        add(colorCheckbox);

        ModifierDensity modifierDensity = new ModifierDensity();
        DeviceSlider densitySlider = new DeviceSlider(sim.getController(), modifierDensity);
        densitySlider.setUnit(dUnit);
        densitySlider.setLabel("Reduced Density");
        densitySlider.setShowBorder(true);
        densitySlider.setMinimum(0);
        densitySlider.setPrecision(3);
        densitySlider.setMaximum(1);
        densitySlider.setNMajor(5);
        densitySlider.setValue(0.04);
        densitySlider.setShowValues(true);
        densitySlider.setEditValues(true);
        add(densitySlider);

        sim.activityIntegrate.setSleepPeriod(0);

        // Simulation Time
        final DisplayTextBox displayCycles = new DisplayTextBox();

        final DataSourceCountTime meterCycles = new DataSourceCountTime(sim.integrator);
        displayCycles.setPrecision(6);
        DataPumpListener pump = new DataPumpListener(meterCycles, displayCycles);
        sim.integrator.getEventManager().addListener(pump);
        displayCycles.setLabel("Simulation time");

        //temperature selector
        tempSlider = new DeviceThermoSlider(sim.getController(), sim.integrator);
        tempSlider.setUnit(tUnit);
//        tempSlider.setPrecision(1);
        tempSlider.setMinimum(0.0);
        tempSlider.setMaximum(600.0);
        tempSlider.setSliderMajorValues(3);
        tempSlider.setAdiabatic();

        add(tempSlider);

        //meters and displays

        DataSourceCountTime timeCounter = new DataSourceCountTime(sim.integrator);

        //add meter and display for current kinetic temperature

        MeterTemperature thermometer = new MeterTemperature(sim.box, space.D());
        DataFork temperatureFork = new DataFork();
        final DataPumpListener temperaturePump = new DataPumpListener(thermometer, temperatureFork);
        sim.integrator.getEventManager().addListener(temperaturePump);
        final AccumulatorAverageCollapsing temperatureAverage = new AccumulatorAverageCollapsing();
        temperatureAverage.setPushInterval(20);
        final AccumulatorHistory temperatureHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
        temperatureHistory.setTimeDataSource(timeCounter);
        temperatureFork.setDataSinks(new IDataSink[]{temperatureAverage, temperatureHistory});
        final DisplayTextBoxesCAE tBox = new DisplayTextBoxesCAE();
        tBox.setAccumulator(temperatureAverage);
        dataStreamPumps.add(temperaturePump);
        tBox.setUnit(tUnit);
        tBox.setLabel("Measured Temperature (K)");
        tBox.setLabelPosition(CompassDirection.NORTH);

        MeterPotentialEnergyFromIntegrator peMeter = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        final AccumulatorHistory peHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
        peHistory.setTimeDataSource(timeCounter);
        final AccumulatorAverageCollapsing peAccumulator = new AccumulatorAverageCollapsing();
        peAccumulator.setPushInterval(2);
        DataFork peFork = new DataFork(new IDataSink[]{peHistory, peAccumulator});
        final DataSinkExcludeOverlap peExcludeOverlap = new DataSinkExcludeOverlap();
        peExcludeOverlap.setDataSink(peFork);
        peExcludeOverlap.numAtoms = sim.box.getLeafList().size();
        final DataPumpListener pePump = new DataPumpListener(peMeter, peExcludeOverlap);
        sim.integrator.getEventManager().addListener(pePump);
        dataStreamPumps.add(pePump);

        // we do this for the scaling by numAtoms rather than for the overlap exclusion
        int numAtoms = sim.box.getLeafList().size();
        pePump.setInterval(numAtoms > 120 ? 1 : 120 / numAtoms);

        final DisplayPlotXChart ePlot = new DisplayPlotXChart();
        peHistory.setDataSink(ePlot.getDataSet().makeDataSink());
        ePlot.setDoLegend(false);

        ePlot.getPlot().setTitle("Energy History (K)");
        ePlot.setLabel("Energy");
        ePlot.setUnit(eUnit);
        ePlot.setXLabel("Simulation Time (ps)");

        final DisplayPlotXChart tPlot = new DisplayPlotXChart();
        temperatureHistory.setDataSink(tPlot.getDataSet().makeDataSink());
        tPlot.setDoLegend(false);

        tPlot.getPlot().setTitle("Temperature History (K)");
        tPlot.setLabel("Temperature");
        tPlot.setUnit(tUnit);
        tPlot.setXLabel("Simulation Time (ps)");

        MeterPressureHard pMeter = new MeterPressureHard(sim.integrator);
        final AccumulatorAverageCollapsing pAccumulator = new AccumulatorAverageCollapsing();
        final DataPumpListener pPump = new DataPumpListener(pMeter, pAccumulator);
        sim.integrator.getEventManager().addListener(pPump);
        pAccumulator.setPushInterval(50);
        dataStreamPumps.add(pPump);

        final DisplayTextBoxesCAE pDisplay = new DisplayTextBoxesCAE();
        pDisplay.setLabel("Reduced Pressure");
        pDisplay.setAccumulator(pAccumulator);
        pDisplay.setUnit(pUnit);
        final DisplayTextBoxesCAE peDisplay = new DisplayTextBoxesCAE();
        peDisplay.setAccumulator(peAccumulator);
        peDisplay.setLabel("Potential Energy (K)");
        peDisplay.setUnit(eUnit);

        meterLargestCluster = new MeterLargestCluster(sim.box);
        DataFork forkCluster = new DataFork();
        DataPumpListener pumpCluster = new DataPumpListener(meterLargestCluster, forkCluster, 10);
        sim.integrator.getEventManager().addListener(pumpCluster);
        AccumulatorHistory clusterHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
        clusterHistory.setTimeDataSource(timeCounter);
        forkCluster.addDataSink(clusterHistory);
        DisplayPlotXChart clusterHistoryPlot = new DisplayPlotXChart();
        clusterHistory.addDataSink(clusterHistoryPlot.getDataSet().makeDataSink());
        clusterHistoryPlot.setLabel("Cluster History");
        clusterHistoryPlot.setDoLegend(false);
        clusterHistoryPlot.setXLabel("Simulation Time (ps)");
        add(clusterHistoryPlot);

        meterClusterSizes = new MeterClusterSizes(sim.box);
        DisplayPlotXChart clusterHistogramPlot = new DisplayPlotXChart();
        DataPumpListener pumpClusterHistogram = new DataPumpListener(meterClusterSizes, clusterHistogramPlot.getDataSet().makeDataSink(), 100);
        sim.integrator.getEventManager().addListener(pumpClusterHistogram);
        clusterHistogramPlot.setLabel("Cluster Histogram");
        clusterHistogramPlot.setDoLegend(false);
        clusterHistogramPlot.getPlot().setXLog(true);
        clusterHistogramPlot.getPlot().setXRange(0, Math.log(sim.box.getLeafList().size()) / Math.log(10));
        clusterHistogramPlot.getPlot().setYRange(0, 1);
        add(clusterHistogramPlot);

        sim.integrator.setTimeStep(0.1);
        sim.integrator.setThermostat(IntegratorMD.ThermostatType.ANDERSEN_SCALING);
        sim.integrator.setThermostatInterval(100);

        //************* Lay out components ****************//

        getDisplayBox(sim.box).setScale(0.7);


        final IAction temperatureAction = new IAction() {
            public void actionPerformed() {
                resetDataAction.actionPerformed();
            }
        };
        tempSlider.setSliderPostAction(temperatureAction);
        tempSlider.setRadioGroupPostAction(temperatureAction);

        final IAction resetAction = new IAction() {
            public void actionPerformed() {
                sim.integrator.reset();

                temperaturePump.actionPerformed();
                tBox.putData(temperatureAverage.getData());

                pMeter.reset();
                pPump.actionPerformed();
                pDisplay.putData(pAccumulator.getData());
                peDisplay.putData(peAccumulator.getData());

                getDisplayBox(sim.box).graphic().repaint();

                displayCycles.putData(meterCycles.getData());
            }
        };

        this.getController().getReinitButton().setPostAction(resetAction);
        this.getController().getResetAveragesButton().setPostAction(resetAction);

        DeviceDelaySlider delaySlider = new DeviceDelaySlider(sim.getController(), sim.activityIntegrate);
        GridBagConstraints vertGBC = SimulationPanel.getVertGBC();

        getPanel().controlPanel.add(delaySlider.graphic(), vertGBC);

        JPanel etPanel = new JPanel(new MigLayout("flowy"));
        etPanel.add(ePlot.getPanel(), "");
        etPanel.add(tPlot.getPanel(), "");
        addAsTab(etPanel, "Energy", true);

        add(displayCycles);
        add(tBox);
        add(pDisplay);
        add(peDisplay);

    }

    protected class ModifierDensity implements Modifier {

        public void setValue(double density) {
            if (density == 0) {
                throw new IllegalArgumentException("density must be positive");
            }
            int N = sim.box.getLeafList().size();
            double rho = N / sim.box.getBoundary().volume();
            double d = Math.pow(density / rho, 1.0 / sim.getSpace().D());
//            if (d > 1) {3
            //assume one type of atom
            ((DiameterHashByType) getDisplayBox(sim.box).getDiameterHash()).setDiameter(sim.species.getLeafType(), d);
            ((P2SquareWell) sim.potentialWrapper.getWrappedPotential()).setCoreDiameter(d);
            new BoxImposePbc(sim.box, space).actionPerformed();
            try {
                sim.integrator.reset();
            } catch (ConfigurationOverlapException e) {
                // can happen when increasing diameter
            }
            meterClusterSizes.clusterer.setNbrMax(d * 1.5);
            meterLargestCluster.clusterer.setNbrMax(d * 1.5);
            getDisplayBox(sim.box).repaint();
        }

        public double getValue() {
            double sigma = ((P2SquareWell) sim.potentialWrapper.getWrappedPotential()).getCoreDiameter();
            int N = sim.box.getLeafList().size();
            double rho = N / sim.box.getBoundary().volume();
            double density = rho * sigma * sigma * sigma;
            return density;
        }

        public Dimension getDimension() {
            return Length.DIMENSION;
        }

        public String getLabel() {
            return "Density";
        }

        public String toString() {
            return getLabel();
        }
    }

    public static class DataSinkExcludeOverlap extends DataProcessor {

        public DataSinkExcludeOverlap() {
            myData = new DataDouble();
        }

        public IData processData(IData data) {
            if (Double.isInfinite(data.getValue(0))) {
                return null;
            }
            myData.E(data);
            myData.TE(1.0 / numAtoms);
            return myData;
        }

        protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
            return inputDataInfo;
        }

        public int numAtoms;
        protected final DataDouble myData;
    }

    public static void main(String[] args) {
        Space space = Space2D.getInstance();
        if (args.length != 0) {
            try {
                int D = Integer.parseInt(args[0]);
                if (D == 3) {
                    space = Space3D.getInstance();
                }
            } catch (NumberFormatException e) {
            }
        }

        NucleationGraphic swmdGraphic = new NucleationGraphic(new Swmd(space));
        SimulationGraphic.makeAndDisplayFrame
                (swmdGraphic.getPanel(), APP_NAME);
    }
}
