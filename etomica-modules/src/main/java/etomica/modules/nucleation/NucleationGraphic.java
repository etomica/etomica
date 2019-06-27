/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.nucleation;

import etomica.action.BoxImposePbc;
import etomica.action.IAction;
import etomica.action.SimulationRestart;
import etomica.atom.DiameterHashByType;
import etomica.config.ConfigurationLattice;
import etomica.data.*;
import etomica.data.histogram.HistogramDiscrete;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterPressureHard;
import etomica.data.meter.MeterTemperature;
import etomica.data.types.DataDouble;
import etomica.exception.ConfigurationOverlapException;
import etomica.graphics.*;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMD;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.modifier.Modifier;
import etomica.modules.swmd.Swmd;
import etomica.potential.P2SquareWell;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space2d.Space2D;
import etomica.space3d.Space3D;
import etomica.units.*;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Length;
import etomica.units.systems.MKS;
import etomica.util.Constants.CompassDirection;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.awt.*;
import java.util.ArrayList;

public class NucleationGraphic extends SimulationGraphic {

    private final static String APP_NAME = "Nucleation";
    private final static int REPAINT_INTERVAL = 1;
    protected DeviceThermoSlider tempSlider;
    protected Unit tUnit, eUnit, dUnit, pUnit, mUnit;
    protected Swmd sim;

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

        eUnit = new UnitRatio(Joule.UNIT, Mole.UNIT);
        mUnit = new UnitRatio(Gram.UNIT, Mole.UNIT);
        if (sim.getSpace().D() == 2) {
            dUnit = new UnitRatio(Mole.UNIT,
                    new MKS().area());
            Unit[] units = new Unit[]{Bar.UNIT, new PrefixedUnit(Prefix.NANO, Meter.UNIT)};
            double[] exponents = new double[]{1.0, 1.0};
            pUnit = new CompoundUnit(units, exponents);
        } else {
            dUnit = new UnitRatio(Mole.UNIT, Liter.UNIT);
            pUnit = Bar.UNIT;
        }

        ((P2SquareWell) sim.potentialWrapper.getWrappedPotential()).setEpsilon(eUnit.toSim(5000));
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
        ModifierDensity modifierDensity = new ModifierDensity();
        DeviceSlider densitySlider = new DeviceSlider(sim.getController(), modifierDensity);
        densitySlider.setLabel("Density");
        densitySlider.setShowBorder(true);
        densitySlider.setMinimum(0);
        densitySlider.setPrecision(3);
        densitySlider.setMaximum(0.1);
        densitySlider.setNMajor(5);
        densitySlider.setValue(400.0 / (100.0 * 100.0));
        add(densitySlider);

        sim.activityIntegrate.setSleepPeriod(0);

        // Simulation Time
        final DisplayTextBox displayCycles = new DisplayTextBox();

        final DataSourceCountTime meterCycles = new DataSourceCountTime(sim.integrator);
        displayCycles.setPrecision(6);
        DataPump pump = new DataPump(meterCycles, displayCycles);
        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(pump));
        displayCycles.setLabel("Simulation time");

        //temperature selector
        tempSlider = new DeviceThermoSlider(sim.getController(), sim.integrator);
        tempSlider.setUnit(tUnit);
//        tempSlider.setPrecision(1);
        tempSlider.setMinimum(0.0);
        tempSlider.setMaximum(1500.0);
        tempSlider.setSliderMajorValues(3);
        tempSlider.setAdiabatic();

        add(tempSlider);

        //display of box, timer
        ColorSchemeByType colorScheme = new ColorSchemeByType();
        colorScheme.setColor(sim.species.getLeafType(), Color.red);
        getDisplayBox(sim.box).setColorScheme(new ColorSchemeByType());

        //meters and displays

        DataSourceCountTime timeCounter = new DataSourceCountTime(sim.integrator);

        //add meter and display for current kinetic temperature

        MeterTemperature thermometer = new MeterTemperature(sim.box, space.D());
        DataFork temperatureFork = new DataFork();
        final DataPump temperaturePump = new DataPump(thermometer, temperatureFork);
        IntegratorListenerAction temperaturePumpListener = new IntegratorListenerAction(temperaturePump);
        sim.integrator.getEventManager().addListener(temperaturePumpListener);
        temperaturePumpListener.setInterval(1);
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
        final DataPump pePump = new DataPump(peMeter, peExcludeOverlap);
        IntegratorListenerAction pePumpListener = new IntegratorListenerAction(pePump);
        sim.integrator.getEventManager().addListener(pePumpListener);
        dataStreamPumps.add(pePump);

        // we do this for the scaling by numAtoms rather than for the overlap exclusion
        int numAtoms = sim.box.getLeafList().size();
        pePumpListener.setInterval(numAtoms > 120 ? 1 : 120 / numAtoms);

        final DisplayPlot ePlot = new DisplayPlot();
        peHistory.setDataSink(ePlot.getDataSet().makeDataSink());
        ePlot.setDoLegend(false);

        ePlot.getPlot().setTitle("Energy History (J/mol)");
        ePlot.setDoLegend(true);
        ePlot.setLabel("Energy");
        ePlot.setUnit(eUnit);
        ePlot.setXUnit(Picosecond.UNIT);

        MeterPressureHard pMeter = new MeterPressureHard(sim.integrator);
        final AccumulatorAverageCollapsing pAccumulator = new AccumulatorAverageCollapsing();
        final DataPumpListener pPump = new DataPumpListener(pMeter, pAccumulator);
        sim.integrator.getEventManager().addListener(pPump);
        pAccumulator.setPushInterval(50);
        dataStreamPumps.add(pPump);

        final DisplayTextBoxesCAE pDisplay = new DisplayTextBoxesCAE();
        pDisplay.setLabel(sim.getSpace().D() == 3 ? "Pressure (bar)" : "Pressure (bar-nm)");
        pDisplay.setAccumulator(pAccumulator);
        pDisplay.setUnit(pUnit);
        final DisplayTextBoxesCAE peDisplay = new DisplayTextBoxesCAE();
        peDisplay.setAccumulator(peAccumulator);
        peDisplay.setLabel("Potential Energy (J/mol)");
        peDisplay.setUnit(eUnit);

        MeterLargestCluster meterLargestCluster = new MeterLargestCluster(sim.box);
        DataFork forkCluster = new DataFork();
        DataPumpListener pumpCluster = new DataPumpListener(meterLargestCluster, forkCluster, 10);
        sim.integrator.getEventManager().addListener(pumpCluster);
        AccumulatorHistory clusterHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
        forkCluster.addDataSink(clusterHistory);
        DisplayPlot clusterHistoryPlot = new DisplayPlot();
        clusterHistory.addDataSink(clusterHistoryPlot.getDataSet().makeDataSink());
        clusterHistoryPlot.setLabel("cluster history");
        add(clusterHistoryPlot);

        AccumulatorHistogram clusterHistogram = new AccumulatorHistogram(new HistogramDiscrete(1e-10));
        forkCluster.addDataSink(clusterHistogram);
        DisplayPlot clusterHistogramPlot = new DisplayPlot();
        clusterHistogram.addDataSink(clusterHistogramPlot.getDataSet().makeDataSink());
        clusterHistogramPlot.setLabel("cluster histogram");
        add(clusterHistogramPlot);

        final DeviceNSelector nSlider = new DeviceNSelector(sim.getController());
        nSlider.setResetAction(new SimulationRestart(sim));
        nSlider.setSpecies(sim.species);
        nSlider.setBox(sim.box);
        nSlider.setMinimum(0);
        nSlider.setMaximum(sim.getSpace().D() == 3 ? 500 : 168);
        nSlider.setLabel("Number of Atoms");
        nSlider.setShowBorder(true);
        nSlider.setShowValues(true);
        // add a listener to adjust the thermostat interval for different
        // system sizes (since we're using ANDERSEN_SINGLE.  Smaller systems
        // don't need as much thermostating.
        ChangeListener nListener = new ChangeListener() {
            public void stateChanged(ChangeEvent evt) {
                final int n = (int) nSlider.getValue() > 0 ? (int) nSlider.getValue() : 1;
                sim.integrator.setThermostatInterval(n > 40 ? 1 : 40 / n);
                peExcludeOverlap.numAtoms = n;

                getDisplayBox(sim.box).repaint();
            }
        };
        nSlider.getSlider().addChangeListener(nListener);
        nListener.stateChanged(null);
        JPanel nSliderPanel = new JPanel(new GridLayout(0, 1));
        nSliderPanel.setBorder(new TitledBorder(null, "Number of Molecules", TitledBorder.CENTER, TitledBorder.TOP));
        nSlider.setShowBorder(false);
        nSlider.setNMajor(4);
        nSliderPanel.add(nSlider.graphic());

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

                // Reset temperature (THIS IS NOT WORKING)
                temperaturePump.actionPerformed();
                tBox.putData(temperatureAverage.getData());
                tBox.repaint();

                // IS THIS WORKING?
                pPump.actionPerformed();
                pDisplay.putData(pAccumulator.getData());
                pDisplay.repaint();
                peDisplay.putData(peAccumulator.getData());
                peDisplay.repaint();

                getDisplayBox(sim.box).graphic().repaint();

                displayCycles.putData(meterCycles.getData());
                displayCycles.repaint();
            }
        };

        this.getController().getReinitButton().setPostAction(new IAction() {
            public void actionPerformed() {
                sim.integrator.setThermostat(IntegratorMD.ThermostatType.ANDERSEN);
                sim.integrator.doThermostat();
                sim.integrator.setThermostat(IntegratorMD.ThermostatType.ANDERSEN_SINGLE);
                resetAction.actionPerformed();
            }
        });
        this.getController().getResetAveragesButton().setPostAction(resetAction);

        DeviceDelaySlider delaySlider = new DeviceDelaySlider(sim.getController(), sim.activityIntegrate);
        GridBagConstraints vertGBC = SimulationPanel.getVertGBC();

        getPanel().controlPanel.add(delaySlider.graphic(), vertGBC);

        add(ePlot);
        add(displayCycles);
        add(tBox);
        add(pDisplay);
        add(peDisplay);

        java.awt.Dimension d = ePlot.getPlot().getPreferredSize();
        d.width -= 50;
        ePlot.getPlot().setSize(d);
    }

    protected class ModifierDensity implements Modifier {

        public void setValue(double density) {
            if (density == 0) throw new IllegalArgumentException("density must be positive");
            int N = sim.box.getLeafList().size();
            double rho = N / sim.box.getBoundary().volume();
            double d = Math.pow(density / rho, 1.0 / sim.getSpace().D());
            if (d > 10.0) {
                throw new IllegalArgumentException("diameter can't exceed 10.0A");
            }
            //assume one type of atom
            ((DiameterHashByType) getDisplayBox(sim.box).getDiameterHash()).setDiameter(sim.species.getLeafType(), d);
            ((P2SquareWell) sim.potentialWrapper.getWrappedPotential()).setCoreDiameter(d);
            new BoxImposePbc(sim.box, space).actionPerformed();
            try {
                sim.integrator.reset();
            } catch (ConfigurationOverlapException e) {
                // can happen when increasing diameter
            }
            getDisplayBox(sim.box).repaint();
        }

        public double getValue() {
            return ((P2SquareWell) sim.potentialWrapper.getWrappedPotential()).getCoreDiameter();
        }

        public Dimension getDimension() {
            return Length.DIMENSION;
        }

        public String getLabel() {
            return "Atom Diameter";
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
