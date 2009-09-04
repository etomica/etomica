package etomica.modules.mu;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemListener;
import java.util.ArrayList;

import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.border.TitledBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import etomica.action.BoxImposePbc;
import etomica.action.IAction;
import etomica.api.IAtomPositioned;
import etomica.api.IAtomTypeSphere;
import etomica.api.IBox;
import etomica.api.IFunction;
import etomica.api.IMolecule;
import etomica.api.IVectorMutable;
import etomica.box.RandomPositionSourceRectangular;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.AccumulatorHistogram;
import etomica.data.AccumulatorHistory;
import etomica.data.DataDump;
import etomica.data.DataFork;
import etomica.data.DataPipe;
import etomica.data.DataProcessor;
import etomica.data.DataProcessorFunction;
import etomica.data.DataPump;
import etomica.data.DataPumpListener;
import etomica.data.DataSourceCountTime;
import etomica.data.DataSplitter;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.meter.MeterDensity;
import etomica.data.meter.MeterNMolecules;
import etomica.data.meter.MeterProfileByVolume;
import etomica.data.meter.MeterWidomInsertion;
import etomica.data.types.DataDouble;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.exception.ConfigurationOverlapException;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DeviceBox;
import etomica.graphics.DeviceDelaySlider;
import etomica.graphics.DeviceNSelector;
import etomica.graphics.DeviceThermoSlider;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTable;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.Drawable;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.integrator.IntegratorBox;
import etomica.modifier.Modifier;
import etomica.modifier.ModifierNMolecule;
import etomica.nbr.list.PotentialMasterList;
import etomica.space.ISpace;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Dimension;
import etomica.units.Energy;
import etomica.units.Fraction;
import etomica.units.Length;
import etomica.units.Null;
import etomica.units.Pixel;
import etomica.util.HistogramDiscrete;
import etomica.util.HistoryCollapsingAverage;

public class MuGraphic extends SimulationGraphic {

    private final static String APP_NAME = "Chemical Potential";
    private final static int REPAINT_INTERVAL = 1;
    protected DeviceThermoSlider tempSlider;
    public ItemListener potentialChooserListener;
    protected Mu sim;

    public MuGraphic(final Mu simulation, ISpace _space) {

    	super(simulation, TABBED_PANE, APP_NAME, REPAINT_INTERVAL, _space, simulation.getController());

        ArrayList<DataPump> dataStreamPumps = getController().getDataStreamPumps();

        final IAction resetDataAction = new IAction() {
            public void actionPerformed() {
                getController().getResetAveragesButton().press();
            }
        };

    	this.sim = simulation;

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
        tempSlider = new DeviceThermoSlider(sim.getController());
        tempSlider.setIsothermalButtonsVisibility(false);
        tempSlider.setPrecision(1);
        tempSlider.setMinimum(0.0);
        tempSlider.setMaximum(10.0);
        tempSlider.setSliderMajorValues(4);
        tempSlider.setAdiabatic();
        tempSlider.setIntegrator(sim.integrator);

        JPanel statePanel = new JPanel(new GridBagLayout());
        GridBagConstraints gbc2 = new GridBagConstraints();
        gbc2.gridx = 0;  gbc2.gridy = 0;
        statePanel.add(tempSlider.graphic(), gbc2);

        GridBagConstraints vertGBC = SimulationPanel.getVertGBC();


        DeviceBox sigABox = new DeviceBox();
        DeviceBox epsABox = new DeviceBox();
        DeviceBox lamABox = new DeviceBox();

        JPanel potentialPanelA = new JPanel(new GridBagLayout());
        JPanel potentialSubPanelA = new JPanel(new GridLayout(0,1));
        potentialSubPanelA.add(sigABox.graphic());
        potentialSubPanelA.add(epsABox.graphic());
        potentialSubPanelA.add(lamABox.graphic());
        potentialPanelA.add(potentialSubPanelA,vertGBC);

        DeviceBox sigBBox = new DeviceBox();
        DeviceBox epsBBox = new DeviceBox();
        DeviceBox lamBBox = new DeviceBox();

        JPanel potentialPanelB = new JPanel(new GridBagLayout());
        JPanel potentialSubPanelB = new JPanel(new GridLayout(0,1));
        potentialSubPanelB.add(sigBBox.graphic());
        potentialSubPanelB.add(epsBBox.graphic());
        potentialSubPanelB.add(lamBBox.graphic());
        potentialPanelB.add(potentialSubPanelB,vertGBC);

        //
        // Tabbed pane for state, potential, controls pages
        //
        JTabbedPane setupPanel = new JTabbedPane();
        setupPanel.add(statePanel, "State");
        setupPanel.add(potentialPanelA, "Potential (A)");
        setupPanel.add(potentialPanelB, "Potential (B)");

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
	    MeterDensity densityMeter = new MeterDensity(sim.getSpace());
        densityMeter.setBox(sim.box);
	    final DisplayTextBox densityBox = new DisplayTextBox();
        final DataPumpListener densityPump = new DataPumpListener(densityMeter, densityBox, 100);
        sim.integrator.getEventManager().addListener(densityPump);
        dataStreamPumps.add(densityPump);
        densityBox.setUnit(Null.UNIT);
	    densityBox.setLabel("Total Density");

        MeterProfileByVolume densityProfileMeterA = new MeterProfileByVolume(space);
        densityProfileMeterA.setBox(sim.box);
        MeterNMolecules meterNMoleculesA = new MeterNMolecules();
        meterNMoleculesA.setSpecies(sim.speciesA);
        densityProfileMeterA.setDataSource(meterNMoleculesA);
        AccumulatorAverageFixed densityProfileAvgA = new AccumulatorAverageFixed(10);
        densityProfileAvgA.setPushInterval(10);
        DataDump profileDumpA = new DataDump();
        densityProfileAvgA.addDataSink(profileDumpA, new AccumulatorAverage.StatType[]{AccumulatorAverage.StatType.AVERAGE});
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
        densityProfileAvgB.addDataSink(profileDumpB, new AccumulatorAverage.StatType[]{AccumulatorAverage.StatType.AVERAGE});
        DataPumpListener profilePumpB = new DataPumpListener(densityProfileMeterB, densityProfileAvgB, 100);
        sim.integrator.getEventManager().addListener(profilePumpB);
        dataStreamPumps.add(profilePumpB);

        DisplayPlot profilePlot = new DisplayPlot();
        densityProfileAvgA.addDataSink(profilePlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{AccumulatorAverage.StatType.AVERAGE});
        densityProfileAvgB.addDataSink(profilePlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{AccumulatorAverage.StatType.AVERAGE});
        profilePlot.setLegend(new DataTag[]{densityProfileMeterA.getTag()}, "A");
        profilePlot.setLegend(new DataTag[]{densityProfileMeterB.getTag()}, "B");
        profilePlot.setLabel("Density");

        MeterDensitySides meterDensityA = new MeterDensitySides(sim.box, sim.speciesA);
        DataSplitter splitterDensityA = new DataSplitter();
        DataPumpListener pumpDensityA = new DataPumpListener(meterDensityA, splitterDensityA);
        sim.integrator.getEventManager().addListener(pumpDensityA);
        dataStreamPumps.add(pumpDensityA);
        AccumulatorAverageCollapsing accumulatorDensityIGA = new AccumulatorAverageCollapsing();
        AccumulatorAverageCollapsing accumulatorDensitySQWA = new AccumulatorAverageCollapsing();
        DataFork densityIGAFork = new DataFork();
        DataFork densitySQWAFork = new DataFork();
        splitterDensityA.setDataSink(0, densityIGAFork);
        densityIGAFork.addDataSink(accumulatorDensityIGA);
        splitterDensityA.setDataSink(1, densitySQWAFork);
        densitySQWAFork.addDataSink(accumulatorDensitySQWA);
        DisplayTextBoxesCAE displayDensityIGA = new DisplayTextBoxesCAE();
        displayDensityIGA.setAccumulator(accumulatorDensityIGA);
        displayDensityIGA.setLabel("IG Phase density (A)");
        displayDensityIGA.setDoShowCurrent(false);
        DisplayTextBoxesCAE displayDensitySQWA = new DisplayTextBoxesCAE();
        displayDensitySQWA.setAccumulator(accumulatorDensitySQWA);
        displayDensitySQWA.setLabel("Real Phase density (A)");
        displayDensitySQWA.setDoShowCurrent(false);
        
        final AccumulatorHistory densityIGAHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
        densityIGAFork.addDataSink(densityIGAHistory);
        densityIGAHistory.setPushInterval(100);
        final AccumulatorHistory densitySQWAHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
        densitySQWAFork.addDataSink(densitySQWAHistory);
        densitySQWAHistory.setPushInterval(100);

        MeterDensitySides meterDensityB = new MeterDensitySides(sim.box, sim.speciesB);
        DataSplitter splitterDensityB = new DataSplitter();
        DataPumpListener pumpDensityB = new DataPumpListener(meterDensityB, splitterDensityB);
        sim.integrator.getEventManager().addListener(pumpDensityB);
        dataStreamPumps.add(pumpDensityB);
        AccumulatorAverageCollapsing accumulatorDensityIGB = new AccumulatorAverageCollapsing();
        AccumulatorAverageCollapsing accumulatorDensitySQWB = new AccumulatorAverageCollapsing();
        DataFork densityIGBFork = new DataFork();
        DataFork densitySQWBFork = new DataFork();
        splitterDensityB.setDataSink(0, densityIGBFork);
        densityIGBFork.addDataSink(accumulatorDensityIGB);
        splitterDensityB.setDataSink(1, densitySQWBFork);
        densitySQWBFork.addDataSink(accumulatorDensitySQWB);
        DisplayTextBoxesCAE displayDensityIGB = new DisplayTextBoxesCAE();
        displayDensityIGB.setAccumulator(accumulatorDensityIGB);
        displayDensityIGB.setLabel("IG Phase density (B)");
        displayDensityIGB.setDoShowCurrent(false);
        DisplayTextBoxesCAE displayDensitySQWB = new DisplayTextBoxesCAE();
        displayDensitySQWB.setAccumulator(accumulatorDensitySQWB);
        displayDensitySQWB.setLabel("Real Phase density (B)");
        displayDensitySQWB.setDoShowCurrent(false);

        final AccumulatorHistory densityIGBHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
        densityIGBFork.addDataSink(densityIGBHistory);
        densityIGBHistory.setPushInterval(100);
        final AccumulatorHistory densitySQWBHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
        densitySQWBFork.addDataSink(densitySQWBHistory);
        densitySQWBHistory.setPushInterval(100);

        MeterWidomInsertion meterMuA = new MeterWidomInsertion(space, sim.getRandom());
        meterMuA.setIntegrator(sim.integrator);
        meterMuA.setNInsert(1);
        meterMuA.setResidual(true);
        meterMuA.setSpecies(sim.speciesA);
        meterMuA.setPositionSource(new RandomPositionSourceRectangular(space, sim.getRandom()) {
            public IVectorMutable randomPosition() {
                IVectorMutable v;
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
        DisplayTextBoxesCAE muDisplayA = new DisplayTextBoxesCAE();
        muDisplayA.setAccumulator(muAvgA);
        muDisplayA.setLabel("exp(-U/kT) (A)");
        muDisplayA.setDoShowCurrent(false);
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
        muHistogramTableA.setColumnHeader(new DataTag[]{((DataInfoFunction)muHistogramA.getDataInfo()).getXDataSource().getIndependentTag()}, "U");
        muHistogramTableA.setColumnHeader(new DataTag[]{muHistogramA.getTag()}, "probability");
        muHistogramTableA.setLabel("Insertion Energy (A)");
        muHistogramTableA.setShowingRowLabels(false);
        
        AccumulatorHistory muHistoryA = new AccumulatorHistory(new HistoryCollapsingAverage());
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
            public IVectorMutable randomPosition() {
                IVectorMutable v;
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
        DisplayTextBoxesCAE muDisplayB = new DisplayTextBoxesCAE();
        muDisplayB.setAccumulator(muAvgB);
        muDisplayB.setLabel("exp(-U/kT) (B)");
        muDisplayB.setDoShowCurrent(false);
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
        muHistogramTableB.setColumnHeader(new DataTag[]{((DataInfoFunction)muHistogramB.getDataInfo()).getXDataSource().getIndependentTag()}, "U");
        muHistogramTableB.setColumnHeader(new DataTag[]{muHistogramB.getTag()}, "probability");
        muHistogramTableB.setLabel("Insertion Energy (B)");
        muHistogramTableB.setShowingRowLabels(false);

        AccumulatorHistory muHistoryB = new AccumulatorHistory(new HistoryCollapsingAverage());
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
        pressureSplitter.setDataSink(0, accumulatorPressureIG);
        pressureSplitter.setDataSink(1, accumulatorPressureSQW);
        DisplayTextBoxesCAE pressureIGDisplay = new DisplayTextBoxesCAE();
        pressureIGDisplay.setAccumulator(accumulatorPressureIG);
        pressureIGDisplay.setLabel("IG Phase Pressure");
        pressureIGDisplay.setDoShowCurrent(false);
        DisplayTextBoxesCAE pressureSQWDisplay = new DisplayTextBoxesCAE();
        pressureSQWDisplay.setAccumulator(accumulatorPressureSQW);
        pressureSQWDisplay.setLabel("Real Phase Pressure");
        pressureSQWDisplay.setDoShowCurrent(false);
        
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
                        IVectorMutable p = ((IAtomPositioned)m.getChildList().getAtom(0)).getPosition();
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
                        IVectorMutable p = ((IAtomPositioned)m.getChildList().getAtom(0)).getPosition();
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


	    ActionListener isothermalListener = new ActionListener() {
	        public void actionPerformed(ActionEvent event) {
                // we can't tell if we're isothermal here...  :(
                // if we're adiabatic, we'll re-set the temperature elsewhere
                resetDataAction.actionPerformed();
            }
        };
		tempSlider.setSliderPostAction(resetDataAction);
        tempSlider.addRadioGroupActionListener(isothermalListener);

        IAction resetAction = new IAction() {
        	public void actionPerformed() {
        	    sim.integrator.reset();

        	    // Reset density (Density is set and won't change, but
        		// do this anyway)
        		densityPump.actionPerformed();
        		densityBox.repaint();

        		getDisplayBox(sim.box).graphic().repaint();
        		
        		displayCycles.putData(meterCycles.getData());
        		displayCycles.repaint();
        	}
        };

        this.getController().getReinitButton().setPostAction(resetAction);
        this.getController().getResetAveragesButton().setPostAction(resetAction);

        DeviceDelaySlider delaySlider = new DeviceDelaySlider(sim.getController(), sim.activityIntegrate);
        
        getPanel().controlPanel.add(setupPanel, vertGBC);
        getPanel().controlPanel.add(delaySlider.graphic(), vertGBC);

        // cause metrics tab to be added
    	add(displayCycles);
    	getPanel().metricPanel.remove(displayCycles.graphic());
    	// now populate metrics tab manually so everything fits
    	JPanel metricSubPanel = new JPanel(new GridBagLayout());
    	GridBagConstraints gbc = new GridBagConstraints();
    	gbc.insets = new Insets(0, 20, 0, 20);
        gbc.gridy = 0;
    	metricSubPanel.add(displayCycles.graphic(), gbc);
        metricSubPanel.add(densityBox.graphic(), gbc);
        gbc.gridy = 1;
        metricSubPanel.add(muDisplayA.graphic(), gbc);
        metricSubPanel.add(displayDensityIGA.graphic(), gbc);
        gbc.gridy = 2;
        metricSubPanel.add(muDisplayB.graphic(), gbc);
        metricSubPanel.add(displayDensitySQWA.graphic(), gbc);
        gbc.gridy = 3;
        metricSubPanel.add(pressureIGDisplay.graphic(), gbc);
        metricSubPanel.add(displayDensityIGB.graphic(), gbc);
        gbc.gridy = 4;
        metricSubPanel.add(pressureSQWDisplay.graphic(), gbc);
        metricSubPanel.add(displayDensitySQWB.graphic(), gbc);
        getPanel().metricPanel.add(metricSubPanel);
        
        add(profilePlot);
        add(muPlot);
        add(muHistogramTableA);
        add(muHistogramTableB);
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

        protected IEtomicaDataInfo processDataInfo(
                IEtomicaDataInfo inputDataInfo) {
            dataInfo = new DataInfoFunction("chemical potential", Energy.DIMENSION, ((DataInfoFunction)inputDataInfo).getXDataSource());
            data = new DataFunction(new int[]{inputDataInfo.getLength()});
            dataInfo.addTag(tag);
            return dataInfo;
        }

        public DataPipe getDataCaster(IEtomicaDataInfo inputDataInfo) {
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
            ((IAtomTypeSphere)species.getAtomType(0)).setDiameter(d);
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

        public DataSinkExcludeOverlap(IBox box) {
            myData = new DataDouble();
            this.box = box;
        }
        
        public DataPipe getDataCaster(IEtomicaDataInfo incomingDataInfo) {
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

        protected IEtomicaDataInfo processDataInfo(IEtomicaDataInfo inputDataInfo) {
            return inputDataInfo;
        }
        
        protected final IBox box;
        protected final DataDouble myData;
    }

    public static void main(String[] args) {
        int dim = 2;
        if (args.length > 0) {
            dim = Integer.parseInt(args[0]);
        }
        ISpace space = Space.getInstance(dim);
        
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
