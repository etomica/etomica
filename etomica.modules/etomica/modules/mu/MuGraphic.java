package etomica.modules.mu;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
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
import etomica.data.DataProcessorChemicalPotential;
import etomica.data.DataProcessorFunction;
import etomica.data.DataPump;
import etomica.data.DataPumpListener;
import etomica.data.DataSourceCountTime;
import etomica.data.DataSourcePositionedBoltzmannFactor;
import etomica.data.DataSplitter;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataSink;
import etomica.data.IEtomicaDataInfo;
import etomica.data.meter.MeterDensity;
import etomica.data.meter.MeterEnergy;
import etomica.data.meter.MeterKineticEnergyFromIntegrator;
import etomica.data.meter.MeterNMolecules;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterProfile;
import etomica.data.meter.MeterProfileByVolume;
import etomica.data.meter.MeterTemperature;
import etomica.data.meter.MeterWidomInsertion;
import etomica.data.types.DataDouble;
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
import etomica.listener.IntegratorListenerAction;
import etomica.modifier.Modifier;
import etomica.modifier.ModifierNMolecule;
import etomica.nbr.list.PotentialMasterList;
import etomica.space.ISpace;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Angstrom;
import etomica.units.Dimension;
import etomica.units.Fraction;
import etomica.units.Length;
import etomica.units.Picosecond;
import etomica.units.Pixel;
import etomica.util.HistogramDiscrete;
import etomica.util.Constants.CompassDirection;

public class MuGraphic extends SimulationGraphic {

    private final static String APP_NAME = "Square-Well Molecular Dynamics";
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
        displayCycles.setLabel("Simulation time");
        
        //temperature selector
        tempSlider = new DeviceThermoSlider(sim.getController());
        tempSlider.setPrecision(1);
        tempSlider.setMinimum(0.0);
        tempSlider.setMaximum(4.0);
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
        setupPanel.add(potentialPanelA, "PotentialA");
        setupPanel.add(potentialPanelB, "PotentialB");

        ModifierAtomDiameter sigModifier = new ModifierAtomDiameter(this, sim.speciesA, sim.potentialAA, sim.potentialAB, sim.potentialBB);
        ModifierEpsilon epsModifier = new ModifierEpsilon(sim.potentialAA, sim.potentialAB, sim.potentialBB, sim.integrator);
        ModifierLambda lamModifier = new ModifierLambda(sim.potentialAA, sim.potentialAB, sim.potentialBB, sim.integrator);
        sigABox.setModifier(sigModifier);
        sigABox.setLabel("Core Diameter ("+Angstrom.UNIT.symbol()+")");
        epsABox.setModifier(epsModifier);
        lamABox.setModifier(lamModifier);
        sigABox.setController(sim.getController());
        epsABox.setController(sim.getController());
        lamABox.setController(sim.getController());

        ModifierAtomDiameter sigBModifier = new ModifierAtomDiameter(this, sim.speciesB, sim.potentialBB, sim.potentialAB, sim.potentialAA);
        ModifierEpsilon epsBModifier = new ModifierEpsilon(sim.potentialBB, sim.potentialAB, sim.potentialAA, sim.integrator);
        ModifierLambda lamBModifier = new ModifierLambda(sim.potentialBB, sim.potentialAB, sim.potentialAA, sim.integrator);
        sigBBox.setModifier(sigBModifier);
        sigBBox.setLabel("Core Diameter ("+Angstrom.UNIT.symbol()+")");
        epsBBox.setModifier(epsBModifier);
        lamBBox.setModifier(lamBModifier);
        sigBBox.setController(sim.getController());
        epsBBox.setController(sim.getController());
        lamBBox.setController(sim.getController());

        //display of box, timer
        ColorSchemeByType colorScheme = new ColorSchemeByType(sim);
        colorScheme.setColor(sim.speciesA.getLeafType(),java.awt.Color.red);
        getDisplayBox(sim.box).setColorScheme(new ColorSchemeByType(sim));

	    //meters and displays
        DataSourceCountTime timeCounter = new DataSourceCountTime(sim.integrator);

        //add meter and display for current kinetic temperature

		MeterTemperature thermometer = new MeterTemperature(sim.box, space.D());
        DataFork temperatureFork = new DataFork();
        final DataPumpListener temperaturePump = new DataPumpListener(thermometer,temperatureFork, 100);
        IntegratorListenerAction temperaturePumpListener = new IntegratorListenerAction(temperaturePump);
        sim.integrator.getEventManager().addListener(temperaturePumpListener);
        temperaturePumpListener.setInterval(1);
        final AccumulatorAverageCollapsing temperatureAverage = new AccumulatorAverageCollapsing();
        temperatureAverage.setPushInterval(20);
        final AccumulatorHistory temperatureHistory = new AccumulatorHistory();
        temperatureHistory.setTimeDataSource(timeCounter);
		temperatureFork.setDataSinks(new IDataSink[]{temperatureAverage,temperatureHistory});
        final DisplayTextBoxesCAE tBox = new DisplayTextBoxesCAE();
        tBox.setAccumulator(temperatureAverage);
		dataStreamPumps.add(temperaturePump);
		tBox.setLabel("Measured Temperature");
		tBox.setLabelPosition(CompassDirection.NORTH);

		// Number density box
	    MeterDensity densityMeter = new MeterDensity(sim.getSpace());
        densityMeter.setBox(sim.box);
	    final DisplayTextBox densityBox = new DisplayTextBox();
        final DataPumpListener densityPump = new DataPumpListener(densityMeter, densityBox, 100);
        sim.integrator.getEventManager().addListener(densityPump);
        dataStreamPumps.add(densityPump);
	    densityBox.setLabel("Density");
	    
		MeterEnergy eMeter = new MeterEnergy(sim.integrator.getPotentialMaster(), sim.box);
        final AccumulatorHistory energyHistory = new AccumulatorHistory();
        energyHistory.setTimeDataSource(timeCounter);
        final DataSinkExcludeOverlap eExcludeOverlap = new DataSinkExcludeOverlap(sim.box);
        eExcludeOverlap.setDataSink(energyHistory);
        final DataPumpListener energyPump = new DataPumpListener(eMeter, eExcludeOverlap, 100);
        sim.integrator.getEventManager().addListener(energyPump);
        dataStreamPumps.add(energyPump);
		
		MeterPotentialEnergyFromIntegrator peMeter = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        final AccumulatorHistory peHistory = new AccumulatorHistory();
        peHistory.setTimeDataSource(timeCounter);
        final AccumulatorAverageCollapsing peAccumulator = new AccumulatorAverageCollapsing();
        peAccumulator.setPushInterval(2);
        DataFork peFork = new DataFork(new IDataSink[]{peHistory, peAccumulator});
        final DataSinkExcludeOverlap peExcludeOverlap = new DataSinkExcludeOverlap(sim.box);
        peExcludeOverlap.setDataSink(peFork);
        final DataPumpListener pePump = new DataPumpListener(peMeter, peExcludeOverlap, 100);
        sim.integrator.getEventManager().addListener(pePump);
        dataStreamPumps.add(pePump);

		MeterKineticEnergyFromIntegrator keMeter = new MeterKineticEnergyFromIntegrator(sim.integrator);
        final AccumulatorHistory keHistory = new AccumulatorHistory();
        keHistory.setTimeDataSource(timeCounter);
        // we do this for the scaling by numAtoms rather than for the overlap exclusion
        final DataSinkExcludeOverlap keExcludeOverlap = new DataSinkExcludeOverlap(sim.box);
        keExcludeOverlap.setDataSink(keHistory);
        final DataPumpListener kePump = new DataPumpListener(keMeter, keExcludeOverlap, 100);
        sim.integrator.getEventManager().addListener(kePump);
        dataStreamPumps.add(kePump);
        
        final DisplayPlot ePlot = new DisplayPlot();
        energyHistory.setDataSink(ePlot.getDataSet().makeDataSink());
        ePlot.setLegend(new DataTag[]{energyHistory.getTag()}, "Total");
        peHistory.setDataSink(ePlot.getDataSet().makeDataSink());
        ePlot.setLegend(new DataTag[]{peHistory.getTag()}, "Potential");
        keHistory.setDataSink(ePlot.getDataSet().makeDataSink());
        ePlot.setLegend(new DataTag[]{keHistory.getTag()}, "Kinetic");

        ePlot.getPlot().setTitle("Energy History (J/mol)");
		ePlot.setDoLegend(true);
		ePlot.setLabel("Energy");
		ePlot.setXUnit(Picosecond.UNIT);
		
//        MeterPressureHard pMeter = new MeterPressureHard(sim.getSpace());
//        pMeter.setIntegrator(sim.integrator);
//        final AccumulatorAverageCollapsing pAccumulator = new AccumulatorAverageCollapsing();
//        final DataPump pPump = new DataPump(pMeter, pAccumulator);
//        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(pPump));
//        pAccumulator.setPushInterval(50);
//        dataStreamPumps.add(pPump);

//        final DisplayTextBoxesCAE pDisplay = new DisplayTextBoxesCAE();
//        pDisplay.setLabel(sim.getSpace().D() == 3 ? "Pressure (bar)" : "Pressure (bar-nm)");
//        pDisplay.setAccumulator(pAccumulator);
        final DisplayTextBoxesCAE peDisplay = new DisplayTextBoxesCAE();
        peDisplay.setAccumulator(peAccumulator);
        peDisplay.setLabel("Potential Energy (J/mol)");
        
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
        profilePlot.setLabel("Density");

        DisplayPlot muPlot = new DisplayPlot();
        MeterProfile muProfileMeter = new MeterProfile(space, sim.getRandom());
        muProfileMeter.setBox(sim.box);
        DataSourcePositionedBoltzmannFactor meterChemicalPotential = new DataSourcePositionedBoltzmannFactor(space);
        meterChemicalPotential.setIntegrator(sim.integrator);
        meterChemicalPotential.setSpecies(sim.speciesA);
        muProfileMeter.setDataSource(meterChemicalPotential);
        AccumulatorAverageFixed chemicalPotentialAverage = new AccumulatorAverageFixed(10);
        chemicalPotentialAverage.setPushInterval(10);
        DataPumpListener muProfilePump = new DataPumpListener(muProfileMeter, chemicalPotentialAverage, 100);
        DataProcessorChemicalPotential dataProcessorChemicalPotential = new DataProcessorChemicalPotential();
        dataProcessorChemicalPotential.setDensityProfileDump(profileDumpA);
        dataProcessorChemicalPotential.setIntegrator(sim.integrator);
        chemicalPotentialAverage.addDataSink(dataProcessorChemicalPotential, new AccumulatorAverage.StatType[]{AccumulatorAverage.StatType.AVERAGE});
        dataProcessorChemicalPotential.setDataSink(muPlot.getDataSet().makeDataSink());
        muPlot.setLegend(new DataTag[]{dataProcessorChemicalPotential.getTag()}, "mu");

        muPlot.setLabel("Chemical Potential");
        muPlot.setDoLegend(false);
        add(muPlot);
        sim.integrator.getEventManager().addListener(muProfilePump);
        dataStreamPumps.add(muProfilePump);
        
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

        
        DataSourceWallPressureMu meterPressure = new DataSourceWallPressureMu(space, sim.p1BoundaryA);
        meterPressure.setIntegrator(sim.integrator);
        DataSplitter pressureSplitter = new DataSplitter();
        DataPumpListener pressurePump = new DataPumpListener(meterPressure, pressureSplitter);
        sim.integrator.getEventManager().addListener(pressurePump);
        AccumulatorAverageCollapsing accumulatorPressureIG = new AccumulatorAverageCollapsing();
        AccumulatorAverageCollapsing accumulatorPressureSQW = new AccumulatorAverageCollapsing();
        pressureSplitter.setDataSink(0, accumulatorPressureIG);
        pressureSplitter.setDataSink(1, accumulatorPressureSQW);
        DisplayTextBoxesCAE pressureIGDisplay = new DisplayTextBoxesCAE();
        pressureIGDisplay.setAccumulator(accumulatorPressureIG);
        pressureIGDisplay.setLabel("Ideal Gas Pressure");
        DisplayTextBoxesCAE pressureSQWDisplay = new DisplayTextBoxesCAE();
        pressureSQWDisplay.setAccumulator(accumulatorPressureSQW);
        pressureSQWDisplay.setLabel("Square-Well Pressure");


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
        nSlider.setMaximum(2000);
        nSlider.setLabel("Number of Atoms");
        nSlider.setShowBorder(true);
        nSlider.setShowValues(true);
        nSlider.setEditValues(true);
        // add a listener to adjust the thermostat interval for different
        // system sizes (since we're using ANDERSEN_SINGLE.  Smaller systems 
        // don't need as much thermostating.
        ChangeListener nListener = new ChangeListener() {
            public void stateChanged(ChangeEvent evt) {
                final int n = (int)nSlider.getValue() > 0 ? (int)nSlider.getValue() : 1;
                sim.integrator.setThermostatInterval(n > 40 ? 1 : 40/n);

                getDisplayBox(sim.box).repaint();
            }
        };
        nSlider.getSlider().addChangeListener(nListener);
        nListener.stateChanged(null);
        JPanel nSliderPanel = new JPanel(new GridLayout(0,1));
        nSliderPanel.setBorder(new TitledBorder(null, "Number of Molecules", TitledBorder.CENTER, TitledBorder.TOP));
        nSlider.setShowBorder(false);
        nSlider.setNMajor(4);
        nSliderPanel.add(nSlider.graphic());
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

        		// Reset temperature (THIS IS NOT WORKING)
                temperaturePump.actionPerformed();
                tBox.putData(temperatureAverage.getData());
                tBox.repaint();

                // IS THIS WORKING?
                peDisplay.putData(peAccumulator.getData());
                peDisplay.repaint();

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

    	add(displayCycles);
    	add(densityBox);
    	add(tBox);
        add(muDisplayA);
        add(muDisplayB);
    	add(peDisplay);
        add(ePlot);
        add(profilePlot);
        add(muHistogramTableA);
        add(muHistogramTableB);
        add(pressureIGDisplay);
        add(pressureSQWDisplay);
    	
        java.awt.Dimension d = ePlot.getPlot().getPreferredSize();
        d.width -= 50;
        ePlot.getPlot().setSize(d);
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
            return "Square-well extent";
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
            if (newValue > 1.75) {
                // our potential neighbor range is 4, so cap lambda at 1.75 (sigma<=2)
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
            return p2.getLambda();
        }

        public Dimension getDimension() {
            return Fraction.DIMENSION;
        }
        
        public String getLabel() {
            return "Square-well extent";
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
            return "Atom Diameter";
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
