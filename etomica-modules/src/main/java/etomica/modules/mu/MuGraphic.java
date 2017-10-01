/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.mu;

import etomica.action.BoxImposePbc;
import etomica.action.IAction;
import etomica.atom.DiameterHashByType;
import etomica.box.Box;
import etomica.box.RandomPositionSourceRectangular;
import etomica.data.*;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.histogram.HistogramDiscrete;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.MeterNMolecules;
import etomica.data.meter.MeterProfileByVolume;
import etomica.data.meter.MeterWidomInsertion;
import etomica.data.types.DataDouble;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.exception.ConfigurationOverlapException;
import etomica.graphics.*;
import etomica.integrator.IntegratorBox;
import etomica.math.function.IFunction;
import etomica.modifier.Modifier;
import etomica.modifier.ModifierNMolecule;
import etomica.molecule.IMolecule;
import etomica.nbr.list.PotentialMasterList;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.SpeciesSpheresMono;
import etomica.units.dimensions.Dimension;
import etomica.units.*;
import etomica.units.dimensions.Energy;
import etomica.units.dimensions.Fraction;
import etomica.units.dimensions.Length;
import etomica.units.dimensions.Null;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.awt.*;
import java.awt.event.ItemListener;
import java.util.ArrayList;

public class MuGraphic extends SimulationGraphic {

    private final static String APP_NAME = "Chemical Potential";
    private final static int REPAINT_INTERVAL = 1;
    protected DeviceThermoSlider tempSlider;
    public ItemListener potentialChooserListener;
    protected Mu sim;

    public MuGraphic(final Mu simulation, Space _space) {

    	super(simulation, TABBED_PANE, APP_NAME, REPAINT_INTERVAL, _space, simulation.getController());

        this.sim = simulation;

        sim.p1BoundaryA.setLongWall(0, false, false);
        sim.p1BoundaryA.setLongWall(0, true, false);
        sim.p1BoundaryA.setBox(sim.box);
        sim.p1BoundaryA.setDrawingThickness(3);
    	getDisplayBox(sim.box).addDrawable(sim.p1BoundaryA);
    	getDisplayBox(sim.box).repaint();
    	getDisplayBox(sim.box).addDrawable(new Drawable() {
            public void draw(Graphics g, int[] origin, double toPixels) {
                int width = (int)(sim.box.getBoundary().getBoxSize().getX(0)*toPixels);
                g.setFont(new Font(null, Font.BOLD, 12));
                g.drawString("Ideal Gas Phase", origin[0]+width/4-50, origin[1]-5);
                g.drawString("Real Phase", origin[0]+3*width/4-35, origin[1]-5);
            }
        });

        ArrayList<DataPump> dataStreamPumps = getController().getDataStreamPumps();

        final IAction resetDataAction = new IAction() {
            public void actionPerformed() {
                getController().getResetAveragesButton().press();
            }
        };

    	getController().getSimRestart().setConfiguration(sim.configuration);

    	getDisplayBox(sim.box).addDrawable(new Drawable() {
            public void draw(Graphics g, int[] origin, double toPixels) {
                int x1 = origin[0]+(int)(0.5*toPixels*sim.box.getBoundary().getBoxSize().getX(0));
                int y1 = origin[1];
                int h = (int)(toPixels*sim.box.getBoundary().getBoxSize().getX(1));
                int w = 2;
                g.setColor(Color.green);
                g.fillRect(x1-w, y1, w, h);
            }
    	});

        getDisplayBox(sim.box).setPixelUnit(new Pixel(40/sim.box.getBoundary().getBoxSize().getX(1)));

        //combo box to select potentials

        // Simulation Time
        final DisplayTextBox displayCycles = new DisplayTextBox();

        final DataSourceCountTime meterCycles = new DataSourceCountTime(sim.integrator);
        displayCycles.setPrecision(6);
        DataPumpListener pump = new DataPumpListener(meterCycles,displayCycles);
        sim.integrator.getEventManager().addListener(pump);
        displayCycles.setUnit(Null.UNIT);
        displayCycles.setLabel("Simulation time");
        
        //temperature selector
        tempSlider = new DeviceThermoSlider(sim.getController(), sim.integrator);
        tempSlider.setIsothermalButtonsVisibility(false);
        tempSlider.setPrecision(1);
        tempSlider.setMinimum(0.0);
        tempSlider.setMaximum(10.0);
        tempSlider.setSliderMajorValues(4);
        tempSlider.setAdiabatic();

        JPanel statePanel = new JPanel(new GridBagLayout());
        GridBagConstraints gbc2 = new GridBagConstraints();
        gbc2.gridx = 0;  gbc2.gridy = 0;
        statePanel.add(tempSlider.graphic(), gbc2);

        GridBagConstraints vertGBC = SimulationPanel.getVertGBC();


        DeviceBox sigABox = new DeviceBox();
        DeviceBox epsABox = new DeviceBox();
        DeviceBox lamABox = new DeviceBox();


        DeviceBox sigBBox = new DeviceBox();
        DeviceBox epsBBox = new DeviceBox();
        DeviceBox lamBBox = new DeviceBox();

        JPanel potentialPanel = new JPanel(new GridBagLayout());

        JPanel potentialPanelA = new JPanel(new GridLayout(0,1));
        potentialPanelA.add(sigABox.graphic());
        potentialPanelA.add(lamABox.graphic());
        potentialPanelA.add(epsABox.graphic());
        TitledBorder border = new TitledBorder("A");
        border.setTitleJustification(TitledBorder.LEFT);
        potentialPanelA.setBorder(border);
        gbc2.gridx = 0; gbc2.gridy = 0;
        potentialPanel.add(potentialPanelA, gbc2);

        JPanel potentialPanelB = new JPanel(new GridLayout(0,1));
        potentialPanelB.add(sigBBox.graphic());
        potentialPanelB.add(lamBBox.graphic());
        potentialPanelB.add(epsBBox.graphic());
        border = new TitledBorder("B");
        border.setTitleJustification(TitledBorder.LEFT);
        potentialPanelB.setBorder(border);
        gbc2.gridx = 1;
        potentialPanel.add(potentialPanelB, gbc2);

        //
        // Tabbed pane for state, potential, controls pages
        //
        JTabbedPane setupPanel = new JTabbedPane();
        setupPanel.add(statePanel, "State");
        setupPanel.add(potentialPanel, "Potential");

        ModifierAtomDiameter sigModifier = new ModifierAtomDiameter(this, sim.speciesA, sim.potentialAA, sim.potentialAB, sim.potentialBB);
        ModifierEpsilon epsModifier = new ModifierEpsilon(sim.potentialAA, sim.potentialAB, sim.potentialBB, sim.integrator);
        ModifierLambda lamModifier = new ModifierLambda(sim.potentialAA, sim.potentialAB, sim.potentialBB, sim.integrator);
        sigABox.setModifier(sigModifier);
        sigABox.setLabel("sigma");
        epsABox.setModifier(epsModifier);
        lamABox.setModifier(lamModifier);
        sigABox.setController(sim.getController());
        epsABox.setController(sim.getController());
        lamABox.setController(sim.getController());

        ModifierAtomDiameter sigBModifier = new ModifierAtomDiameter(this, sim.speciesB, sim.potentialBB, sim.potentialAB, sim.potentialAA);
        ModifierEpsilon epsBModifier = new ModifierEpsilon(sim.potentialBB, sim.potentialAB, sim.potentialAA, sim.integrator);
        ModifierLambda lamBModifier = new ModifierLambda(sim.potentialBB, sim.potentialAB, sim.potentialAA, sim.integrator);
        sigBBox.setModifier(sigBModifier);
        sigBBox.setLabel("sigma");
        epsBBox.setModifier(epsBModifier);
        lamBBox.setModifier(lamBModifier);
        sigBBox.setController(sim.getController());
        epsBBox.setController(sim.getController());
        lamBBox.setController(sim.getController());

        //display of box, timer
        ColorSchemeByType colorScheme = new ColorSchemeByType(sim);
        colorScheme.setColor(sim.speciesA.getLeafType(),java.awt.Color.red);
        getDisplayBox(sim.box).setColorScheme(new ColorSchemeByType(sim));

		// Number density box

        MeterProfileByVolume densityProfileMeterA = new MeterProfileByVolume(space);
        densityProfileMeterA.setBox(sim.box);
        MeterNMolecules meterNMoleculesA = new MeterNMolecules();
        meterNMoleculesA.setSpecies(sim.speciesA);
        densityProfileMeterA.setDataSource(meterNMoleculesA);
        AccumulatorAverageFixed densityProfileAvgA = new AccumulatorAverageFixed(10);
        densityProfileAvgA.setPushInterval(10);
        DataDump profileDumpA = new DataDump();
        densityProfileAvgA.addDataSink(profileDumpA, new AccumulatorAverage.StatType[]{densityProfileAvgA.AVERAGE});
        DataPumpListener profilePumpA = new DataPumpListener(densityProfileMeterA, densityProfileAvgA, 100);
        sim.integrator.getEventManager().addListener(profilePumpA);
        dataStreamPumps.add(profilePumpA);

        MeterProfileByVolume densityProfileMeterB = new MeterProfileByVolume(space);
        densityProfileMeterB.setBox(sim.box);
        MeterNMolecules meterNMoleculesB = new MeterNMolecules();
        meterNMoleculesB.setSpecies(sim.speciesB);
        densityProfileMeterB.setDataSource(meterNMoleculesB);
        AccumulatorAverageFixed densityProfileAvgB = new AccumulatorAverageFixed(10);
        densityProfileAvgB.setPushInterval(10);
        DataDump profileDumpB = new DataDump();
        densityProfileAvgB.addDataSink(profileDumpB, new AccumulatorAverage.StatType[]{densityProfileAvgB.AVERAGE});
        DataPumpListener profilePumpB = new DataPumpListener(densityProfileMeterB, densityProfileAvgB, 100);
        sim.integrator.getEventManager().addListener(profilePumpB);
        dataStreamPumps.add(profilePumpB);

        DisplayPlot profilePlot = new DisplayPlot();
        densityProfileAvgA.addDataSink(profilePlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{densityProfileAvgA.AVERAGE});
        densityProfileAvgB.addDataSink(profilePlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{densityProfileAvgB.AVERAGE});
        profilePlot.setLegend(new DataTag[]{densityProfileMeterA.getTag()}, "A");
        profilePlot.setLegend(new DataTag[]{densityProfileMeterB.getTag()}, "B");
        profilePlot.setLabel("Density");

        MeterDensitySides meterDensitySQWA = new MeterDensitySides(sim.box, sim.speciesA, false);
        DataFork densitySQWAFork = new DataFork();
        DataPumpListener pumpDensitySQWA = new DataPumpListener(meterDensitySQWA, densitySQWAFork);
        sim.integrator.getEventManager().addListener(pumpDensitySQWA);
        dataStreamPumps.add(pumpDensitySQWA);
        AccumulatorAverageCollapsing accumulatorDensitySQWA = new AccumulatorAverageCollapsing();
        densitySQWAFork.addDataSink(accumulatorDensitySQWA);

        MeterDensitySides meterDensitySQWB = new MeterDensitySides(sim.box, sim.speciesB, false);
        DataFork densitySQWBFork = new DataFork();
        DataPumpListener pumpDensitySQWB = new DataPumpListener(meterDensitySQWB, densitySQWBFork);
        sim.integrator.getEventManager().addListener(pumpDensitySQWB);
        dataStreamPumps.add(pumpDensitySQWB);
        AccumulatorAverageCollapsing accumulatorDensitySQWB = new AccumulatorAverageCollapsing();
        densitySQWBFork.addDataSink(accumulatorDensitySQWB);

        MeterDensitySides meterDensityIGA = new MeterDensitySides(sim.box, sim.speciesA, true);
        DataFork densityIGAFork = new DataFork();
        DataPumpListener pumpDensityIGA = new DataPumpListener(meterDensityIGA, densityIGAFork);
        sim.integrator.getEventManager().addListener(pumpDensityIGA);
        dataStreamPumps.add(pumpDensityIGA);
        AccumulatorAverageCollapsing accumulatorDensityIGA = new AccumulatorAverageCollapsing();
        densityIGAFork.addDataSink(accumulatorDensityIGA);

        MeterDensitySides meterDensityIGB = new MeterDensitySides(sim.box, sim.speciesB, true);
        DataFork densityIGBFork = new DataFork();
        DataPumpListener pumpDensityIGB = new DataPumpListener(meterDensityIGB, densityIGBFork);
        sim.integrator.getEventManager().addListener(pumpDensityIGB);
        dataStreamPumps.add(pumpDensityIGB);
        AccumulatorAverageCollapsing accumulatorDensityIGB = new AccumulatorAverageCollapsing();
        densityIGBFork.addDataSink(accumulatorDensityIGB);

        final AccumulatorHistory densityIGAHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
        densityIGAHistory.setTimeDataSource(meterCycles);
        densityIGAFork.addDataSink(densityIGAHistory);
        densityIGAHistory.setPushInterval(100);
        final AccumulatorHistory densitySQWAHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
        densitySQWAHistory.setTimeDataSource(meterCycles);
        densitySQWAFork.addDataSink(densitySQWAHistory);
        densitySQWAHistory.setPushInterval(100);

        final AccumulatorHistory densityIGBHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
        densityIGBHistory.setTimeDataSource(meterCycles);
        densityIGBFork.addDataSink(densityIGBHistory);
        densityIGBHistory.setPushInterval(100);
        final AccumulatorHistory densitySQWBHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
        densitySQWBHistory.setTimeDataSource(meterCycles);
        densitySQWBFork.addDataSink(densitySQWBHistory);
        densitySQWBHistory.setPushInterval(100);

        MeterWidomInsertion meterMuA = new MeterWidomInsertion(space, sim.getRandom());
        meterMuA.setIntegrator(sim.integrator);
        meterMuA.setNInsert(1);
        meterMuA.setResidual(true);
        meterMuA.setSpecies(sim.speciesA);
        meterMuA.setPositionSource(new RandomPositionSourceRectangular(space, sim.getRandom()) {
            public Vector randomPosition() {
                Vector v;
                do {
                    v = super.randomPosition();
                }
                while (v.getX(0) < 0);
                return v;
            }
        });
        DataFork muForkA = new DataFork();
        DataPumpListener muPumpA = new DataPumpListener(meterMuA, muForkA);
        AccumulatorAverageCollapsing muAvgA = new AccumulatorAverageCollapsing();
        muForkA.addDataSink(muAvgA);
        sim.integrator.getEventManager().addListener(muPumpA);
        dataStreamPumps.add(muPumpA);
        muAvgA.setPushInterval(100);
        DataProcessor uProcessorA = new DataProcessorFunction(new IFunction() {
            public double f(double x) {
                if (x==0) return Double.POSITIVE_INFINITY;
                return -Math.log(x)*sim.integrator.getTemperature();
            }
        });
        muForkA.addDataSink(uProcessorA);
        AccumulatorHistogram muHistogramA = new AccumulatorHistogram(new HistogramDiscrete(1e-10));
        uProcessorA.setDataSink(muHistogramA);
        DisplayTable muHistogramTableA = new DisplayTable();
        muHistogramA.setDataSink(muHistogramTableA.getDataTable().makeDataSink());
        muHistogramTableA.setColumnHeader(new DataTag[]{((DataInfoFunction)muHistogramA.getDataInfo()).getXDataSource().getIndependentTag()}, "E");
        muHistogramTableA.setColumnHeader(new DataTag[]{muHistogramA.getTag()}, "probability");
        muHistogramTableA.setLabel("Insertion Energy (A)");
        muHistogramTableA.setShowingRowLabels(false);
        
        AccumulatorHistory muHistoryA = new AccumulatorHistory(new HistoryCollapsingAverage());
        muHistoryA.setTimeDataSource(meterCycles);
        muForkA.addDataSink(muHistoryA);
        DataProcessorMu muSQWA = new DataProcessorMu(muHistoryA, sim.integrator);
        densitySQWAHistory.setDataSink(muSQWA);
        DataProcessorMu muIGA = new DataProcessorMu(null, sim.integrator);
        densityIGAHistory.setDataSink(muIGA);

        MeterWidomInsertion meterMuB = new MeterWidomInsertion(space, sim.getRandom());
        meterMuB.setIntegrator(sim.integrator);
        meterMuB.setNInsert(1);
        meterMuB.setResidual(true);
        meterMuB.setSpecies(sim.speciesB);
        meterMuB.setPositionSource(new RandomPositionSourceRectangular(space, sim.getRandom()) {
            public Vector randomPosition() {
                Vector v;
                do {
                    v = super.randomPosition();
                }
                while (v.getX(0) < 0);
                return v;
            }
        });
        DataFork muForkB = new DataFork();
        DataPumpListener muPumpB = new DataPumpListener(meterMuB, muForkB);
        AccumulatorAverageCollapsing muAvgB = new AccumulatorAverageCollapsing();
        muForkB.addDataSink(muAvgB);
        sim.integrator.getEventManager().addListener(muPumpB);
        dataStreamPumps.add(muPumpB);
        muAvgB.setPushInterval(100);
        DataProcessor uProcessorB = new DataProcessorFunction(new IFunction() {
            public double f(double x) {
                if (x==0) return Double.POSITIVE_INFINITY;
                return -Math.log(x)*sim.integrator.getTemperature();
            }
        });
        muForkB.addDataSink(uProcessorB);
        AccumulatorHistogram muHistogramB = new AccumulatorHistogram(new HistogramDiscrete(1e-10));
        uProcessorB.setDataSink(muHistogramB);
        DisplayTable muHistogramTableB = new DisplayTable();
        muHistogramB.setDataSink(muHistogramTableB.getDataTable().makeDataSink());
        muHistogramTableB.setColumnHeader(new DataTag[]{((DataInfoFunction)muHistogramB.getDataInfo()).getXDataSource().getIndependentTag()}, "E");
        muHistogramTableB.setColumnHeader(new DataTag[]{muHistogramB.getTag()}, "probability");
        muHistogramTableB.setLabel("Insertion Energy (B)");
        muHistogramTableB.setShowingRowLabels(false);

        AccumulatorHistory muHistoryB = new AccumulatorHistory(new HistoryCollapsingAverage());
        muHistoryB.setTimeDataSource(meterCycles);
        muForkB.addDataSink(muHistoryB);
        DataProcessorMu muSQWB = new DataProcessorMu(muHistoryB, sim.integrator);
        densitySQWBHistory.setDataSink(muSQWB);
        DataProcessorMu muIGB = new DataProcessorMu(null, sim.integrator);
        densityIGBHistory.setDataSink(muIGB);

        DisplayPlot muPlot = new DisplayPlot();
        muPlot.setLabel("Chemical Potential");
        muSQWA.setDataSink(muPlot.getDataSet().makeDataSink());
        muHistoryA.setPushInterval(100);
        muPlot.setLegend(new DataTag[]{muSQWA.getTag()}, "Real Phase (A)");
        muSQWB.setDataSink(muPlot.getDataSet().makeDataSink());
        muHistoryB.setPushInterval(100);
        muPlot.setLegend(new DataTag[]{muSQWB.getTag()}, "Real Phase (B)");
        muIGA.setDataSink(muPlot.getDataSet().makeDataSink());
        muPlot.setLegend(new DataTag[]{muIGA.getTag()}, "IG Phase (A)");
        muIGB.setDataSink(muPlot.getDataSet().makeDataSink());
        muPlot.setLegend(new DataTag[]{muIGB.getTag()}, "IG Phase (B)");

        DataSourceWallPressureMu meterPressure = new DataSourceWallPressureMu(space);
        meterPressure.setIntegrator(sim.integrator);
        DataSplitter pressureSplitter = new DataSplitter();
        DataPumpListener pressurePump = new DataPumpListener(meterPressure, pressureSplitter);
        sim.integrator.getEventManager().addListener(pressurePump);
        dataStreamPumps.add(pressurePump);
        AccumulatorAverageCollapsing accumulatorPressureIG = new AccumulatorAverageCollapsing();
        AccumulatorAverageCollapsing accumulatorPressureSQW = new AccumulatorAverageCollapsing();
        pressureSplitter.setDataSink(0, accumulatorPressureSQW);
        pressureSplitter.setDataSink(1, accumulatorPressureIG);
        
        DisplayTable metricsTable = new DisplayTable();
        metricsTable.setTransposed(true);
        muAvgA.addDataSink(metricsTable.getDataTable().makeDataSink(), new StatType[]{muAvgA.AVERAGE, muAvgA.ERROR});
        muAvgB.addDataSink(metricsTable.getDataTable().makeDataSink(), new StatType[]{muAvgB.AVERAGE, muAvgB.ERROR});
        accumulatorDensityIGA.addDataSink(metricsTable.getDataTable().makeDataSink(), new StatType[]{accumulatorDensityIGA.AVERAGE, accumulatorDensityIGA.ERROR});
        accumulatorDensityIGB.addDataSink(metricsTable.getDataTable().makeDataSink(), new StatType[]{accumulatorDensityIGB.AVERAGE, accumulatorDensityIGB.ERROR});
        accumulatorDensitySQWA.addDataSink(metricsTable.getDataTable().makeDataSink(), new StatType[]{accumulatorDensitySQWA.AVERAGE, accumulatorDensitySQWA.ERROR});
        accumulatorDensitySQWB.addDataSink(metricsTable.getDataTable().makeDataSink(), new StatType[]{accumulatorDensitySQWB.AVERAGE, accumulatorDensitySQWB.ERROR});
        accumulatorPressureIG.addDataSink(metricsTable.getDataTable().makeDataSink(), new StatType[]{accumulatorPressureIG.AVERAGE, accumulatorPressureIG.ERROR});
        accumulatorPressureSQW.addDataSink(metricsTable.getDataTable().makeDataSink(), new StatType[]{accumulatorPressureSQW.AVERAGE, accumulatorPressureSQW.ERROR});
        metricsTable.setColumnHeader(new DataTag[]{muAvgA.getTag()}, "exp(-E/kT) (A)");
        metricsTable.setColumnHeader(new DataTag[]{muAvgB.getTag()}, "exp(-E/kT) (B)");
        metricsTable.setColumnHeader(new DataTag[]{accumulatorDensityIGA.getTag()}, "IG Phase Density (A)");
        metricsTable.setColumnHeader(new DataTag[]{accumulatorDensityIGB.getTag()}, "IG Phase Density (B)");
        metricsTable.setColumnHeader(new DataTag[]{accumulatorDensitySQWA.getTag()}, "Real Phase Density (A)");
        metricsTable.setColumnHeader(new DataTag[]{accumulatorDensitySQWB.getTag()}, "Real Phase Density (B)");
        metricsTable.setColumnHeader(new DataTag[]{accumulatorPressureIG.getTag()}, "IG Phase Pressure");
        metricsTable.setColumnHeader(new DataTag[]{accumulatorPressureSQW.getTag()}, "Real Phase Pressure");
        metricsTable.setRowLabels(new String[]{"Average", "Error"});
        metricsTable.setShowingRowLabels(true);
        metricsTable.setLabel("Metrics");
        
        final DeviceNSelector nSlider = new DeviceNSelector(sim.getController());
        nSlider.setSpecies(sim.speciesA);
        nSlider.setBox(sim.box);
        nSlider.setModifier(new ModifierNMolecule(sim.box, sim.speciesA) {
            public void setValue(double newValue) {
                int d = (int)newValue;
                int oldValue = box.getNMolecules(species);
                if (d < oldValue) {
                    box.setNMolecules(species, d);
                }
                else {
                    for (int i=0; i<(d-oldValue); i++) {
                        IMolecule m = species.makeMolecule();
                        Vector p = m.getChildList().getAtom(0).getPosition();
                        p.setX(0, -7.5);
                        box.addMolecule(m);
                    }
                }
                sim.integrator.reset();
            }
        });
        nSlider.setMinimum(0);
        nSlider.setMaximum(1000);
        nSlider.setNMajor(4);
        nSlider.setShowBorder(true);
        nSlider.setLabel("(A)");
        nSlider.setShowValues(true);
        nSlider.setEditValues(true);
        ChangeListener nListener = new ChangeListener() {
            public void stateChanged(ChangeEvent evt) {
                getDisplayBox(sim.box).repaint();
            }
        };
        nSlider.getSlider().addChangeListener(nListener);

        final DeviceNSelector nSliderB = new DeviceNSelector(sim.getController());
        nSliderB.setSpecies(sim.speciesB);
        nSliderB.setBox(sim.box);
        nSliderB.setModifier(new ModifierNMolecule(sim.box, sim.speciesB) {
            public void setValue(double newValue) {
                int d = (int)newValue;
                int oldValue = box.getNMolecules(species);
                if (d < oldValue) {
                    box.setNMolecules(species, d);
                }
                else {
                    for (int i=0; i<(d-oldValue); i++) {
                        IMolecule m = species.makeMolecule();
                        Vector p = m.getChildList().getAtom(0).getPosition();
                        p.setX(0, -7.5);
                        box.addMolecule(m);
                    }
                }
                sim.integrator.reset();
            }
        });
        nSliderB.setMinimum(0);
        nSliderB.setMaximum(1000);
        nSliderB.setNMajor(4);
        nSliderB.setShowBorder(true);
        nSliderB.setLabel("(B)");
        nSliderB.setShowValues(true);
        nSliderB.setEditValues(true);
        nSliderB.getSlider().addChangeListener(nListener);
        
        JPanel nSliderPanel = new JPanel(new GridLayout(0,1));
        nSliderPanel.setBorder(new TitledBorder(null, "Number of Molecules", TitledBorder.CENTER, TitledBorder.TOP));
        nSliderPanel.add(nSlider.graphic());
        nSliderPanel.add(nSliderB.graphic());
        gbc2.gridx = 0;  gbc2.gridy = 1;
        statePanel.add(nSliderPanel, gbc2);

        //************* Lay out components ****************//

        getDisplayBox(sim.box).setScale(0.7);


		tempSlider.setSliderPostAction(resetDataAction);
        tempSlider.setRadioGroupPostAction(resetDataAction);

        IAction resetAction = new IAction() {
        	public void actionPerformed() {
        	    sim.integrator.reset();

        		getDisplayBox(sim.box).graphic().repaint();
        		
        		displayCycles.putData(meterCycles.getData());
        		displayCycles.repaint();
        	}
        };

        this.getController().getReinitButton().setPostAction(resetAction);
        this.getController().getResetAveragesButton().setPostAction(resetAction);

        DeviceDelaySlider delaySlider = new DeviceDelaySlider(sim.getController(), sim.activityIntegrate);
        
        getPanel().controlPanel.add(displayCycles.graphic(), vertGBC);
        getPanel().controlPanel.add(setupPanel, vertGBC);
        getPanel().controlPanel.add(delaySlider.graphic(), vertGBC);

        add(profilePlot);
        add(muPlot);
        add(muHistogramTableA);
        add(muHistogramTableB);
        add(metricsTable);
    }

    public static class DataProcessorMu extends DataProcessor {
        private static final long serialVersionUID = 1L;
        protected final AccumulatorHistory muHistory;
        protected DataFunction data;
        protected DataInfoFunction dataInfo;
        protected final IntegratorBox integrator;

        public DataProcessorMu(AccumulatorHistory muHistory, IntegratorBox integrator) {
            this.muHistory = muHistory;
            this.integrator = integrator;
        }

        protected IData processData(IData inputData) {
            double[] x = data.getData();
            double temp = integrator.getTemperature();
            IData muData = muHistory == null ? null : muHistory.getData();
            for (int i=0; i<inputData.getLength(); i++) {
                double density = inputData.getValue(i);
                double aexp = muData == null ? 1 : muData.getValue(i);
                if (density*aexp == 0) {
                    x[i] = Double.NaN;
                }
                else {
                    x[i] = temp * Math.log(density/aexp);
                }
            }
            return data;
        }

        protected IDataInfo processDataInfo(
                IDataInfo inputDataInfo) {
            dataInfo = new DataInfoFunction("chemical potential", Energy.DIMENSION, ((DataInfoFunction)inputDataInfo).getXDataSource());
            data = new DataFunction(new int[]{inputDataInfo.getLength()});
            dataInfo.addTag(tag);
            return dataInfo;
        }

        public DataPipe getDataCaster(IDataInfo inputDataInfo) {
            return null;
        }
    }

    public static class ModifierLambda implements Modifier {
        public ModifierLambda(P2SquareWellOneSide p2, P2SquareWellOneSide p2Mix,
                P2SquareWellOneSide p2Other, IntegratorBox integrator) {
            this.integrator = integrator;
            this.p2 = p2;
            this.p2Mix = p2Mix;
            this.p2Other = p2Other;
        }

        public void setValue(double newValue) {
            if (newValue > 1.75) {
                // our potential neighbor range is 4, so cap lambda at 1.75 (sigma<=2)
                throw new IllegalArgumentException();
            }
            p2.setLambda(newValue);
            double sigma = p2.getCoreDiameter();
            double otherLambda = p2Other.getLambda();
            double otherSigma = p2Other.getCoreDiameter();
            p2Mix.setLambda((sigma*newValue+otherSigma*otherLambda)/(sigma+otherSigma));
            ((PotentialMasterList)integrator.getPotentialMaster()).reset();
            try {
                integrator.reset();
            }
            catch (ConfigurationOverlapException e){
                // could already be overlapped from increasing diameter
            }
        }

        public double getValue() {
            return p2.getLambda();
        }

        public Dimension getDimension() {
            return Fraction.DIMENSION;
        }
        
        public String getLabel() {
            return "lambda";
        }


        private static final long serialVersionUID = 1L;
        protected final IntegratorBox integrator;
        protected final P2SquareWellOneSide p2, p2Mix, p2Other;
    }

    public static class ModifierEpsilon implements Modifier {
        public ModifierEpsilon(P2SquareWellOneSide p2, P2SquareWellOneSide p2Mix,
                P2SquareWellOneSide p2Other, IntegratorBox integrator) {
            this.integrator = integrator;
            this.p2 = p2;
            this.p2Mix = p2Mix;
            this.p2Other = p2Other;
        }

        public void setValue(double newValue) {
            if (newValue > 10 || newValue  < 0) {
                throw new IllegalArgumentException();
            }
            p2.setEpsilon(newValue);
            double otherEpsilon = p2Other.getEpsilon();
            p2Mix.setEpsilon(Math.sqrt(newValue*otherEpsilon));
            try {
                integrator.reset();
            }
            catch (ConfigurationOverlapException e){
                // could already be overlapped from increasing diameter
            }
        }

        public double getValue() {
            return p2.getEpsilon();
        }

        public Dimension getDimension() {
            return Fraction.DIMENSION;
        }
        
        public String getLabel() {
            return "epsilon";
        }


        private static final long serialVersionUID = 1L;
        protected final IntegratorBox integrator;
        protected final P2SquareWellOneSide p2, p2Mix, p2Other;
    }

    protected static class ModifierAtomDiameter implements Modifier {

        public ModifierAtomDiameter(MuGraphic simGraphic, SpeciesSpheresMono species, P2SquareWellOneSide p2,
                P2SquareWellOneSide p2Mix, P2SquareWellOneSide p2Other) {
            this.simGraphic = simGraphic;
            this.species = species;
            this.p2 = p2;
            this.p2Mix = p2Mix;
            this.p2Other = p2Other;
        }

        public void setValue(double d) {
            if (d > 2.0) {
                throw new IllegalArgumentException("diameter can't exceed 2.0A");
            }
            //assume one type of atom
            Mu sim = simGraphic.sim;
            ((DiameterHashByType)simGraphic.getDisplayBox(sim.box).getDiameterHash()).setDiameter(species.getLeafType(), d);
            p2.setCoreDiameter(d);
            double sigmaOther = p2Other.getCoreDiameter();
            p2Mix.setCoreDiameter(0.5*(sigmaOther+d));
            new BoxImposePbc(sim.box, sim.getSpace()).actionPerformed();
            ((PotentialMasterList)sim.integrator.getPotentialMaster()).reset();
            try {
                sim.integrator.reset();
            }
            catch (ConfigurationOverlapException e){
                // can happen when increasing diameter
            }
            simGraphic.getDisplayBox(sim.box).repaint();
        }

        public double getValue() {
            return p2.getCoreDiameter();
        }

        public Dimension getDimension() {
            return Length.DIMENSION;
        }
        
        public String getLabel() {
            return "sigma";
        }

        protected final MuGraphic simGraphic;
        protected final SpeciesSpheresMono species;
        protected final P2SquareWellOneSide p2, p2Mix, p2Other;
    }
    
    public static class DataSinkExcludeOverlap extends DataProcessor {

        public DataSinkExcludeOverlap(Box box) {
            myData = new DataDouble();
            this.box = box;
        }
        
        public DataPipe getDataCaster(IDataInfo incomingDataInfo) {
            return null;
        }
        
        public IData processData(IData data) {
            if (Double.isInfinite(data.getValue(0))) {
                return null;
            }
            myData.E(data);
            int numAtoms = box.getLeafList().getAtomCount();
            myData.TE(1.0/numAtoms);
            return myData;
        }

        protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
            return inputDataInfo;
        }
        
        protected final Box box;
        protected final DataDouble myData;
    }

    public static void main(String[] args) {
        int dim = 2;
        if (args.length > 0) {
            dim = Integer.parseInt(args[0]);
        }
        Space space = Space.getInstance(dim);
        
        MuGraphic swmdGraphic = new MuGraphic(new Mu(space), space);
		SimulationGraphic.makeAndDisplayFrame
		        (swmdGraphic.getPanel(), APP_NAME);
    }
    
    public static class Applet extends javax.swing.JApplet {

        public void init() {
	        getRootPane().putClientProperty(
	                        "defeatSystemEventQueueCheck", Boolean.TRUE);
            String dimStr = getParameter("dim");
            int dim = 2;
            if (dimStr != null) {
                dim = Integer.valueOf(dimStr).intValue();
            }
            Space sp = Space.getInstance(dim);
            MuGraphic swmdGraphic = new MuGraphic(new Mu(sp), sp);

		    getContentPane().add(swmdGraphic.getPanel());
	    }

        private static final long serialVersionUID = 1L;
    }

}
