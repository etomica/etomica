package etomica.modules.swmd;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.event.ItemListener;
import java.util.ArrayList;

import javax.swing.JComboBox;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.border.TitledBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import etomica.action.Action;
import etomica.action.BoxImposePbc;
import etomica.action.SimulationRestart;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.AtomTypeSphere;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.AccumulatorHistory;
import etomica.data.Data;
import etomica.data.DataFork;
import etomica.data.DataPipe;
import etomica.data.DataProcessor;
import etomica.data.DataPump;
import etomica.data.DataSink;
import etomica.data.DataSourceCountTime;
import etomica.data.DataSourceFunction;
import etomica.data.DataSourceRmsVelocity;
import etomica.data.DataSourceUniform;
import etomica.data.DataTag;
import etomica.data.IDataInfo;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.DataSourceUniform.LimitType;
import etomica.data.meter.MeterDensity;
import etomica.data.meter.MeterEnergy;
import etomica.data.meter.MeterKineticEnergyFromIntegrator;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterPressureHard;
import etomica.data.meter.MeterRDF;
import etomica.data.meter.MeterTemperature;
import etomica.data.types.DataDouble;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.exception.ConfigurationOverlapException;
import etomica.graphics.ActionConfigWindow;
import etomica.graphics.ActionVelocityWindow;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DeviceBox;
import etomica.graphics.DeviceButton;
import etomica.graphics.DeviceDelaySlider;
import etomica.graphics.DeviceNSelector;
import etomica.graphics.DeviceThermoSlider;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.modifier.Modifier;
import etomica.modifier.ModifierGeneral;
import etomica.potential.P2HardSphere;
import etomica.potential.P2Ideal;
import etomica.potential.P2SquareWell;
import etomica.space.Space;
import etomica.space2d.Space2D;
import etomica.space3d.Space3D;
import etomica.statmech.MaxwellBoltzmann;
import etomica.units.Angstrom;
import etomica.units.Bar;
import etomica.units.CompoundUnit;
import etomica.units.Dalton;
import etomica.units.Dimension;
import etomica.units.DimensionRatio;
import etomica.units.Joule;
import etomica.units.Kelvin;
import etomica.units.Length;
import etomica.units.Liter;
import etomica.units.Meter;
import etomica.units.Mole;
import etomica.units.Null;
import etomica.units.Picosecond;
import etomica.units.Pixel;
import etomica.units.Prefix;
import etomica.units.PrefixedUnit;
import etomica.units.Time;
import etomica.units.Unit;
import etomica.units.UnitRatio;
import etomica.units.systems.MKS;
import etomica.util.DoubleRange;
import etomica.util.HistogramSimple;
import etomica.util.Constants.CompassDirection;

public class SwmdGraphic extends SimulationGraphic {

    private final static String APP_NAME = "Square-Well Molecular Dynamics";
    private final static int REPAINT_INTERVAL = 1;
    protected DeviceThermoSlider tempSlider;
    public ItemListener potentialChooserListener;
    public JComboBox potentialChooser;
    protected P2HardSphere potentialHS;
    protected P2SquareWell potentialSW;
    protected P2Ideal potentialIdeal;
    public DeviceBox sigBox, epsBox, lamBox, massBox;
    public double lambda, epsilon, mass, sigma;
    protected Unit eUnit, dUnit, pUnit;
    protected Swmd sim;
    
    private boolean showConfig = true;

    public SwmdGraphic(final Swmd simulation) {

    	super(simulation, TABBED_PANE, APP_NAME, REPAINT_INTERVAL);

        ArrayList dataStreamPumps = getController().getDataStreamPumps();

    	this.sim = simulation;

        Unit tUnit = Kelvin.UNIT;

        eUnit = new UnitRatio(Joule.UNIT, Mole.UNIT);
        lambda = 1.5;
        epsilon = eUnit.toSim(1500.0);
        mass = 40;
        sigma = 4.0;
        if (sim.getSpace().D() == 2) {
            dUnit = new UnitRatio(Mole.UNIT, 
                                    new MKS().area());
            Unit[] units = new Unit[] {Bar.UNIT, new PrefixedUnit(Prefix.NANO, Meter.UNIT)};
            double[] exponents = new double[] {1.0, 1.0};
            pUnit = new CompoundUnit(units, exponents);
        }
        else {
            dUnit = new UnitRatio(Mole.UNIT, Liter.UNIT);
            pUnit = Bar.UNIT;

        }
        

        if (sim.getSpace().D() == 2) {
            getDisplayBox(sim.box).setPixelUnit(new Pixel(400/sim.box.getBoundary().getDimensions().x(1)));
        }
        else {
            getDisplayBox(sim.box).setPixelUnit(new Pixel(40/sim.box.getBoundary().getDimensions().x(1)));
        }

        sim.activityIntegrate.setSleepPeriod(0);
       
        //combo box to select potentials
        final String idealGas = "Ideal gas";
        final String repulsionOnly = "Repulsion only";
        final String repulsionAttraction = "Repulsion and attraction";
        potentialChooser = new javax.swing.JComboBox(new String[] {
                idealGas, repulsionOnly, repulsionAttraction});

        sigBox = new DeviceBox();
        epsBox = new DeviceBox();
        lamBox = new DeviceBox();
        massBox = new DeviceBox();
        // Unselectable because "Ideal gas" is selected initially
        potentialChooser.setSelectedIndex(0);
        epsBox.setEditable(false);
        lamBox.setEditable(false);

        // Simulation Time
        final DisplayTextBox displayCycles = new DisplayTextBox();

        final DataSourceCountTime meterCycles = new DataSourceCountTime(sim.integrator);
        displayCycles.setPrecision(6);
        DataPump pump= new DataPump(meterCycles,displayCycles);
        sim.integrator.addIntervalAction(pump);
        displayCycles.setLabel("Simulation time");
        
        //temperature selector
        tempSlider = new DeviceThermoSlider(sim.getController());
        tempSlider.setUnit(Kelvin.UNIT);
//        tempSlider.setPrecision(1);
        tempSlider.setMinimum(0.0);
        tempSlider.setMaximum(500.0);
        tempSlider.setSliderMajorValues(3);
        tempSlider.setUnit(tUnit);
        tempSlider.setAdiabatic();
        tempSlider.setIntegrator(sim.integrator);

        JPanel statePanel = new JPanel(new GridBagLayout());
        GridBagConstraints gbc2 = new GridBagConstraints();
        gbc2.gridx = 0;  gbc2.gridy = 0;
        statePanel.add(tempSlider.graphic(), gbc2);

        GridBagConstraints vertGBC = SimulationPanel.getVertGBC();

        JPanel potentialPanel = new JPanel(new GridBagLayout());
        potentialPanel.add(potentialChooser,vertGBC);
        potentialPanel.setBorder(new TitledBorder(null, "Potential Selection", TitledBorder.CENTER, TitledBorder.TOP));
        JPanel parameterPanel = new JPanel(new GridLayout(0,1));
        parameterPanel.add(sigBox.graphic());
        parameterPanel.add(epsBox.graphic());
        parameterPanel.add(lamBox.graphic());
        parameterPanel.add(massBox.graphic());
        potentialPanel.add(parameterPanel,vertGBC);
        
        //
        // Tabbed pane for state, potential, controls pages
        //
        JTabbedPane setupPanel = new JTabbedPane();
        setupPanel.add(statePanel, "State");
        setupPanel.add(potentialPanel, "Potential");

        potentialSW = new P2SquareWell(sim.getSpace(),sigma,lambda,epsilon,true);
        potentialHS = new P2HardSphere(sim.getSpace(),sigma,true);
        potentialIdeal = new P2Ideal(sim.getSpace());
        
        if(potentialChooserListener != null) potentialChooser.removeItemListener(potentialChooserListener);
        
        potentialChooserListener = new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent evt) {
                if(evt.getStateChange() == java.awt.event.ItemEvent.DESELECTED) return; 
                setPotential((String)evt.getItem());
                if((String)evt.getItem() == idealGas ||
                   (String)evt.getItem() == repulsionOnly) {
                    epsBox.setEditable(false);
                    lamBox.setEditable(false);
                }
                else {
                    epsBox.setEditable(true);
                    lamBox.setEditable(true);   
                }
            }
        };
        potentialChooser.addItemListener(potentialChooserListener);
        setPotential((String)potentialChooser.getSelectedItem());

        ModifierAtomDiameter sigModifier = new ModifierAtomDiameter();
        sigModifier.setValue(sigma);
        ModifierGeneral epsModifier = new ModifierGeneral(potentialSW, "epsilon");
        ModifierGeneral lamModifier = new ModifierGeneral(potentialSW, "lambda");
        ModifierGeneral massModifier = new ModifierGeneral(((AtomTypeLeaf)sim.species.getMoleculeType()).getElement(),"mass");
        sigBox.setModifier(sigModifier);
        sigBox.setLabel("Core Diameter ("+Angstrom.UNIT.symbol()+")");
        epsBox.setUnit(eUnit);
        epsBox.setModifier(epsModifier);
        lamBox.setModifier(lamModifier);
        massBox.setModifier(massModifier);
        massBox.setUnit(Dalton.UNIT);
        sigBox.setController(sim.getController());
        epsBox.setController(sim.getController());
        lamBox.setController(sim.getController());
        massBox.setController(sim.getController());

        //display of box, timer
        ColorSchemeByType colorScheme = new ColorSchemeByType();
        colorScheme.setColor(sim.species.getMoleculeType(),java.awt.Color.red);
        getDisplayBox(sim.box).setColorScheme(new ColorSchemeByType());

	    //meters and displays
        final MeterRDF rdfMeter = new MeterRDF(sim.getSpace());
        sim.integrator.addIntervalAction(rdfMeter);
        sim.integrator.setActionInterval(rdfMeter, 1);
        rdfMeter.getXDataSource().setXMax(12.0);
        rdfMeter.setBox(sim.box);
        DisplayPlot rdfPlot = new DisplayPlot();
        DataPump rdfPump = new DataPump(rdfMeter,rdfPlot.getDataSet().makeDataSink());
        sim.integrator.addIntervalAction(rdfPump);
        sim.integrator.setActionInterval(rdfPump, 1);
        
        rdfPlot.setDoLegend(false);
        rdfPlot.getPlot().setTitle("Radial Distribution Function");
        rdfPlot.setLabel("RDF");

		//velocity distribution
        double vMin = 0;
        double vMax = 10;
        DataSourceRmsVelocity meterVelocity = new DataSourceRmsVelocity(new HistogramSimple(100,new DoubleRange(0,vMax)));
        meterVelocity.setIterator(new AtomIteratorLeafAtoms(sim.box));
        AccumulatorAverage rmsAverage = new AccumulatorAverageFixed(10);
        DataPump velocityPump = new DataPump(meterVelocity, rmsAverage);
        sim.integrator.addIntervalAction(velocityPump);
        sim.integrator.setActionInterval(velocityPump, 1);
        rmsAverage.setPushInterval(10);
        dataStreamPumps.add(velocityPump);
        
        final DisplayPlot vPlot = new DisplayPlot();
        rmsAverage.addDataSink(vPlot.getDataSet().makeDataSink(), new StatType[]{StatType.AVERAGE});
        vPlot.setDoLegend(false);
        vPlot.getPlot().setTitle("Velocity Distribution");
        vPlot.setDoLegend(true);
        vPlot.setLabel("Velocity");
		
		final MaxwellBoltzmann.Distribution mbDistribution = new MaxwellBoltzmann.Distribution(sim.getSpace(),sim.integrator.getTemperature(),((AtomTypeLeaf)sim.species.getMoleculeType()).getMass());
		final DataSourceFunction mbSource = new DataSourceFunction("Maxwell Boltzmann Distribution",
                Null.DIMENSION, mbDistribution, 100, "Speed (A/ps)", new DimensionRatio(Length.DIMENSION,Time.DIMENSION));
		DataSourceUniform mbX = mbSource.getXSource();
		mbX.setTypeMax(LimitType.HALF_STEP);
		mbX.setTypeMin(LimitType.HALF_STEP);
		mbX.setNValues(((DataInfoFunction)meterVelocity.getDataInfo()).getLength());
		mbX.setXMin(vMin);
		mbX.setXMax(vMax);
		mbSource.update();
        DataPump mbPump = new DataPump(mbSource,vPlot.getDataSet().makeDataSink());
        sim.integrator.addIntervalAction(mbPump);
        sim.integrator.setActionInterval(mbPump, 1);
		
        DataSourceCountTime timeCounter = new DataSourceCountTime(sim.integrator);

        //add meter and display for current kinetic temperature

		MeterTemperature thermometer = new MeterTemperature();
        thermometer.setBox(sim.box);
        DataFork temperatureFork = new DataFork();
        final DataPump temperaturePump = new DataPump(thermometer,temperatureFork);
        sim.integrator.addIntervalAction(temperaturePump);
        sim.integrator.setActionInterval(temperaturePump, 1);
        final AccumulatorAverageCollapsing temperatureAverage = new AccumulatorAverageCollapsing();
        temperatureAverage.setPushInterval(20);
        final AccumulatorHistory temperatureHistory = new AccumulatorHistory();
        temperatureHistory.setTimeDataSource(timeCounter);
		temperatureFork.setDataSinks(new DataSink[]{temperatureAverage,temperatureHistory});
        final DisplayTextBoxesCAE tBox = new DisplayTextBoxesCAE();
        tBox.setAccumulator(temperatureAverage);
		dataStreamPumps.add(temperaturePump);
        tBox.setUnit(tUnit);
		tBox.setLabel("Measured Temperature (K)");
		tBox.setLabelPosition(CompassDirection.NORTH);

		// Number density box
	    MeterDensity densityMeter = new MeterDensity(sim.getSpace());
        densityMeter.setBox(sim.box);
	    final DisplayTextBox densityBox = new DisplayTextBox();
	    densityBox.setUnit(dUnit);
        final DataPump densityPump = new DataPump(densityMeter, densityBox);
        sim.integrator.addIntervalAction(densityPump);
        sim.integrator.setActionInterval(densityPump, 1);
        dataStreamPumps.add(densityPump);
	    densityBox.setLabel("Density");
	    
		MeterEnergy eMeter = new MeterEnergy(sim.integrator.getPotential());
        eMeter.setBox(sim.box);
        final AccumulatorHistory energyHistory = new AccumulatorHistory();
        energyHistory.setTimeDataSource(timeCounter);
        final DataSinkExcludeOverlap eExcludeOverlap = new DataSinkExcludeOverlap();
        eExcludeOverlap.setDataSink(energyHistory);
        final DataPump energyPump = new DataPump(eMeter, eExcludeOverlap);
        sim.integrator.addIntervalAction(energyPump);
        dataStreamPumps.add(energyPump);
		
		MeterPotentialEnergyFromIntegrator peMeter = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        final AccumulatorHistory peHistory = new AccumulatorHistory();
        peHistory.setTimeDataSource(timeCounter);
        final AccumulatorAverageCollapsing peAccumulator = new AccumulatorAverageCollapsing();
        peAccumulator.setPushInterval(2);
        DataFork peFork = new DataFork(new DataSink[]{peHistory, peAccumulator});
        final DataSinkExcludeOverlap peExcludeOverlap = new DataSinkExcludeOverlap();
        peExcludeOverlap.setDataSink(peFork);
        final DataPump pePump = new DataPump(peMeter, peExcludeOverlap);
        sim.integrator.addIntervalAction(pePump);
        dataStreamPumps.add(pePump);

		MeterKineticEnergyFromIntegrator keMeter = new MeterKineticEnergyFromIntegrator(sim.integrator);
        final AccumulatorHistory keHistory = new AccumulatorHistory();
        keHistory.setTimeDataSource(timeCounter);
        // we do this for the scaling by numAtoms rather than for the overlap exclusion
        final DataSinkExcludeOverlap keExcludeOverlap = new DataSinkExcludeOverlap();
        keExcludeOverlap.setDataSink(keHistory);
        final DataPump kePump = new DataPump(keMeter, keExcludeOverlap);
        sim.integrator.addIntervalAction(kePump);
        dataStreamPumps.add(kePump);
        int numAtoms = sim.box.getLeafList().getAtomCount();
        sim.integrator.setActionInterval(energyPump, numAtoms > 120 ? 1 : 120/numAtoms);
        sim.integrator.setActionInterval(kePump, numAtoms > 120 ? 1 : 120/numAtoms);
        sim.integrator.setActionInterval(pePump, numAtoms > 120 ? 1 : 120/numAtoms);
        
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
		
        MeterPressureHard pMeter = new MeterPressureHard(sim.getSpace());
        pMeter.setIntegrator(sim.integrator);
        final AccumulatorAverageCollapsing pAccumulator = new AccumulatorAverageCollapsing();
        final DataPump pPump = new DataPump(pMeter, pAccumulator);
        sim.integrator.addIntervalAction(pPump);
        pAccumulator.setPushInterval(50);
        dataStreamPumps.add(pPump);

        final DisplayTextBoxesCAE pDisplay = new DisplayTextBoxesCAE();
        pDisplay.setLabel(sim.getSpace().D() == 3 ? "Pressure (bar)" : "Pressure (bar-nm)");
        pDisplay.setAccumulator(pAccumulator);
        pDisplay.setUnit(pUnit);
        final DisplayTextBoxesCAE peDisplay = new DisplayTextBoxesCAE();
        peDisplay.setAccumulator(peAccumulator);
        peDisplay.setLabel("Potential Energy (J/mol)");

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
                final int n = (int)nSlider.getValue() > 0 ? (int)nSlider.getValue() : 1;
                sim.integrator.setThermostatInterval(n > 40 ? 1 : 40/n);
                eExcludeOverlap.numAtoms = n;
                peExcludeOverlap.numAtoms = n;
                keExcludeOverlap.numAtoms = n;

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


	    ChangeListener temperatureListener = new ChangeListener() {
		    public void stateChanged(ChangeEvent event) {

		        mbDistribution.setTemperature(Kelvin.UNIT.toSim(tempSlider.getTemperature()));
		        mbSource.update();
		        vPlot.doUpdate();
		        vPlot.repaint();
		    }
		};
		tempSlider.addTemperatureSliderListener(temperatureListener);

        // show config button
        DeviceButton configButton = new DeviceButton(sim.getController());
        configButton.setLabel("Show Config");
        configButton.setAction(new ActionConfigWindow(sim.box));
        DeviceButton velocityButton = new DeviceButton(sim.getController());
        velocityButton.setLabel("Show Velocities");
        velocityButton.setAction(new ActionVelocityWindow(sim.box));

        Action resetAction = new Action() {
        	public void actionPerformed() {
        	    sim.integrator.initialize();

        	    rdfMeter.reset();

        	    // Reset density (Density is set and won't change, but
        		// do this anyway)
        		densityPump.actionPerformed();
        		densityBox.repaint();

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

        this.getController().getReinitButton().setPostAction(resetAction);
        this.getController().getResetAveragesButton().setPostAction(resetAction);

        DeviceDelaySlider delaySlider = new DeviceDelaySlider(sim.getController(), sim.activityIntegrate);
        
        getPanel().controlPanel.add(setupPanel, vertGBC);
        getPanel().controlPanel.add(delaySlider.graphic(), vertGBC);
        if(showConfig == true) {
            add(configButton);
            add(velocityButton);
        }

    	add(rdfPlot);
    	add(vPlot);
    	add(ePlot);
    	add(displayCycles);
    	add(densityBox);
    	add(tBox);
    	add(pDisplay);
    	add(peDisplay);
    }

    public void setPotential(String potentialDesc) {
        final boolean HS = potentialDesc.equals("Repulsion only"); 
        final boolean SW = potentialDesc.equals("Repulsion and attraction"); 
        sim.getController().doActionNow( new Action() {
            public void actionPerformed() {
                if (HS) {
                    potentialHS.setBox(sim.box);
                    sim.potentialWrapper.setPotential(potentialHS);
                }
                else if (SW) {
                    potentialSW.setBox(sim.box);
                    sim.potentialWrapper.setPotential(potentialSW);
                }
                else {
                    potentialIdeal.setBox(sim.box);
                    sim.potentialWrapper.setPotential(potentialIdeal);
                }
                // wrap atoms back inside the box so we can reassign atoms to cells if needed
                if (sim.integrator.isInitialized()) {
                    try {
                        sim.integrator.reset();
                    } catch(ConfigurationOverlapException e) {}
                }
            }
        });
    }

    protected class ModifierAtomDiameter implements Modifier {

        public void setValue(double d) {
            if (d > 4.0) {
                throw new IllegalArgumentException("diameter can't exceed 4.0A");
            }
            //assume one type of atom
            ((AtomTypeSphere)sim.species.getMoleculeType()).setDiameter(d);
            SwmdGraphic.this.potentialHS.setCollisionDiameter(d);
            SwmdGraphic.this.potentialSW.setCoreDiameter(d);
            new BoxImposePbc(sim.box).actionPerformed();
            try {
                sim.integrator.reset();
            }
            catch (ConfigurationOverlapException e){
                // can happen when increasing diameter
            }
            sigma = d;
            getDisplayBox(sim.box).repaint();
        }

        public double getValue() {
            return sigma;
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
        
        public DataPipe getDataCaster(IDataInfo incomingDataInfo) {
            return null;
        }
        
        public Data processData(Data data) {
            if (Double.isInfinite(data.getValue(0))) {
                return null;
            }
            myData.E(data);
            myData.TE(1.0/numAtoms);
            return myData;
        }

        protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
            return inputDataInfo;
        }
        
        public int numAtoms;
        protected final DataDouble myData;
    }

    public static void main(String[] args) {
        Space space = Space3D.getInstance();
        if(args.length != 0) {
            try {
                int D = Integer.parseInt(args[0]);
                if (D == 3) {
                    space = Space3D.getInstance();
                }
            } catch(NumberFormatException e) {}
        }

        SwmdGraphic swmdGraphic = new SwmdGraphic(new Swmd(space));
		SimulationGraphic.makeAndDisplayFrame
		        (swmdGraphic.getPanel(), APP_NAME);
    }
    
    public static class Applet extends javax.swing.JApplet {

        public void init() {
	        getRootPane().putClientProperty(
	                        "defeatSystemEventQueueCheck", Boolean.TRUE);
            SwmdGraphic swmdGraphic = new SwmdGraphic(new Swmd(Space3D.getInstance()));

		    getContentPane().add(swmdGraphic.getPanel());
	    }

        private static final long serialVersionUID = 1L;
    }

}


