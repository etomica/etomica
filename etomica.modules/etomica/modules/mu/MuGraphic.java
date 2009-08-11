package etomica.modules.mu;

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
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.listener.IntegratorListenerAction;
import etomica.modifier.Modifier;
import etomica.modifier.ModifierGeneral;
import etomica.modifier.ModifierNMolecule;
import etomica.nbr.list.PotentialMasterList;
import etomica.space.ISpace;
import etomica.space.Space;
import etomica.units.Angstrom;
import etomica.units.Dimension;
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
    public DeviceBox sigBox, epsBox, lamBox;
    public double lambda, epsilon, sigma;
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

        lambda = sim.potentialSW.getLambda();
        epsilon = sim.potentialSW.getEpsilon();
        sigma = sim.potentialSW.getCoreDiameter();

        getDisplayBox(sim.box).setPixelUnit(new Pixel(40/sim.box.getBoundary().getBoxSize().getX(1)));

        //combo box to select potentials
        sigBox = new DeviceBox();
        epsBox = new DeviceBox();
        lamBox = new DeviceBox();

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

        JPanel potentialPanel = new JPanel(new GridBagLayout());
        JPanel parameterPanel = new JPanel(new GridLayout(0,1));
        parameterPanel.add(sigBox.graphic());
        parameterPanel.add(epsBox.graphic());
        parameterPanel.add(lamBox.graphic());
        potentialPanel.add(parameterPanel,vertGBC);
        
        //
        // Tabbed pane for state, potential, controls pages
        //
        JTabbedPane setupPanel = new JTabbedPane();
        setupPanel.add(statePanel, "State");
        setupPanel.add(potentialPanel, "Potential");

        ModifierAtomDiameter sigModifier = new ModifierAtomDiameter();
        sigModifier.setValue(sigma);
        ModifierGeneral epsModifier = new ModifierGeneral(new Object[]{sim.potentialSW,sim.p1Wall}, "epsilon");
        ModifierGeneral lamModifier = new ModifierGeneral(new Object[]{sim.potentialSW,sim.p1Wall}, "lambda") {
            public void setValue(double newValue) {
                if (sim.potentialSW.getCoreDiameter()*newValue > 3) {
                    // our potential neighbor range is 4, so cap interaction range at 3
                    throw new IllegalArgumentException();
                }
                super.setValue(newValue);
                ((PotentialMasterList)sim.integrator.getPotentialMaster()).reset();
                try {
                    sim.integrator.reset();
                }
                catch (ConfigurationOverlapException e){
                    // could already be overlapped from increasing diameter
                }
            }
        };
        sigBox.setModifier(sigModifier);
        sigBox.setLabel("Core Diameter ("+Angstrom.UNIT.symbol()+")");
        epsBox.setModifier(epsModifier);
        lamBox.setModifier(lamModifier);
        sigBox.setController(sim.getController());
        epsBox.setController(sim.getController());
        lamBox.setController(sim.getController());

        //display of box, timer
        ColorSchemeByType colorScheme = new ColorSchemeByType(sim);
        colorScheme.setColor(sim.species.getLeafType(),java.awt.Color.red);
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
        final DataSinkExcludeOverlap eExcludeOverlap = new DataSinkExcludeOverlap();
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
        final DataSinkExcludeOverlap peExcludeOverlap = new DataSinkExcludeOverlap();
        peExcludeOverlap.setDataSink(peFork);
        final DataPumpListener pePump = new DataPumpListener(peMeter, peExcludeOverlap, 100);
        sim.integrator.getEventManager().addListener(pePump);
        dataStreamPumps.add(pePump);

		MeterKineticEnergyFromIntegrator keMeter = new MeterKineticEnergyFromIntegrator(sim.integrator);
        final AccumulatorHistory keHistory = new AccumulatorHistory();
        keHistory.setTimeDataSource(timeCounter);
        // we do this for the scaling by numAtoms rather than for the overlap exclusion
        final DataSinkExcludeOverlap keExcludeOverlap = new DataSinkExcludeOverlap();
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
        
        MeterProfileByVolume densityProfileMeter = new MeterProfileByVolume(space);
        densityProfileMeter.setBox(sim.box);
        MeterNMolecules meterNMolecules = new MeterNMolecules();
        meterNMolecules.setSpecies(sim.species);
        densityProfileMeter.setDataSource(meterNMolecules);
        AccumulatorAverageFixed densityProfileAvg = new AccumulatorAverageFixed(10);
        densityProfileAvg.setPushInterval(10);
        DataDump profileDump = new DataDump();
        densityProfileAvg.addDataSink(profileDump, new AccumulatorAverage.StatType[]{AccumulatorAverage.StatType.AVERAGE});
        DataPumpListener profilePump = new DataPumpListener(densityProfileMeter, densityProfileAvg, 100);
        sim.integrator.getEventManager().addListener(profilePump);
        dataStreamPumps.add(profilePump);

        DisplayPlot profilePlot = new DisplayPlot();
        densityProfileAvg.addDataSink(profilePlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{AccumulatorAverage.StatType.AVERAGE});
        profilePlot.setDoLegend(false);
        profilePlot.setLabel("Density");

        DisplayPlot muPlot = new DisplayPlot();
        MeterProfile muProfileMeter = new MeterProfile(space, sim.getRandom());
        muProfileMeter.setBox(sim.box);
        DataSourcePositionedBoltzmannFactor meterChemicalPotential = new DataSourcePositionedBoltzmannFactor(space);
        meterChemicalPotential.setIntegrator(sim.integrator);
        meterChemicalPotential.setSpecies(sim.species);
        muProfileMeter.setDataSource(meterChemicalPotential);
        AccumulatorAverageFixed chemicalPotentialAverage = new AccumulatorAverageFixed(10);
        chemicalPotentialAverage.setPushInterval(10);
        DataPumpListener muProfilePump = new DataPumpListener(muProfileMeter, chemicalPotentialAverage, 100);
        DataProcessorChemicalPotential dataProcessorChemicalPotential = new DataProcessorChemicalPotential();
        dataProcessorChemicalPotential.setDensityProfileDump(profileDump);
        dataProcessorChemicalPotential.setIntegrator(sim.integrator);
        chemicalPotentialAverage.addDataSink(dataProcessorChemicalPotential, new AccumulatorAverage.StatType[]{AccumulatorAverage.StatType.AVERAGE});
        dataProcessorChemicalPotential.setDataSink(muPlot.getDataSet().makeDataSink());
        muPlot.setLegend(new DataTag[]{dataProcessorChemicalPotential.getTag()}, "mu");

        muPlot.setLabel("Chemical Potential");
        muPlot.setDoLegend(false);
        add(muPlot);
        sim.integrator.getEventManager().addListener(muProfilePump);
        dataStreamPumps.add(muProfilePump);
        
        MeterWidomInsertion meterMu = new MeterWidomInsertion(space, sim.getRandom());
        meterMu.setIntegrator(sim.integrator);
        meterMu.setNInsert(1);
        meterMu.setResidual(true);
        meterMu.setSpecies(sim.species);
        meterMu.setPositionSource(new RandomPositionSourceRectangular(space, sim.getRandom()) {
            public IVectorMutable randomPosition() {
                IVectorMutable v;
                do {
                    v = super.randomPosition();
                }
                while (v.getX(0) < 0);
                return v;
            }
        });
        DataFork muFork = new DataFork();
        DataPumpListener muPump = new DataPumpListener(meterMu, muFork);
        AccumulatorAverageCollapsing muAvg = new AccumulatorAverageCollapsing();
        muFork.addDataSink(muAvg);
        sim.integrator.getEventManager().addListener(muPump);
        DisplayTextBoxesCAE muDisplay = new DisplayTextBoxesCAE();
        muDisplay.setAccumulator(muAvg);
        muAvg.setPushInterval(100);
        DataProcessor uProcessor = new DataProcessorFunction(new IFunction() {
            public double f(double x) {
                if (x==0) return Double.POSITIVE_INFINITY;
                return -Math.log(x)*sim.integrator.getTemperature();
            }
        });
        muFork.addDataSink(uProcessor);
        AccumulatorHistogram muHistogram = new AccumulatorHistogram(new HistogramDiscrete(1e-10));
        uProcessor.setDataSink(muHistogram);
        DisplayTable muHistogramTable = new DisplayTable();
        muHistogram.setDataSink(muHistogramTable.getDataTable().makeDataSink());
        muHistogramTable.setColumnHeader(new DataTag[]{((DataInfoFunction)muHistogram.getDataInfo()).getXDataSource().getIndependentTag()}, "U");
        muHistogramTable.setColumnHeader(new DataTag[]{muHistogram.getTag()}, "probability");
        muHistogramTable.setLabel("Insertion Energy");

        final DeviceNSelector nSlider = new DeviceNSelector(sim.getController());
        nSlider.setSpecies(sim.species);
        nSlider.setBox(sim.box);
        nSlider.setModifier(new ModifierNMolecule(sim.box, sim.species) {
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
        add(muDisplay);
    	add(peDisplay);
        add(ePlot);
        add(profilePlot);
        add(muHistogramTable);
    	
        java.awt.Dimension d = ePlot.getPlot().getPreferredSize();
        d.width -= 50;
        ePlot.getPlot().setSize(d);
    }

    protected class ModifierAtomDiameter implements Modifier {

        public void setValue(double d) {
            if (d > 3.0 || d*sim.potentialSW.getLambda() > 3.0) {
                throw new IllegalArgumentException("diameter can't exceed 3.0A");
            }
            //assume one type of atom
            ((IAtomTypeSphere)sim.species.getLeafType()).setDiameter(d);
            sim.potentialSW.setCoreDiameter(d);
            sim.p1Wall.setSigma(d);
            new BoxImposePbc(sim.box, space).actionPerformed();
            ((PotentialMasterList)sim.integrator.getPotentialMaster()).reset();
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
        
        public DataPipe getDataCaster(IEtomicaDataInfo incomingDataInfo) {
            return null;
        }
        
        public IData processData(IData data) {
            if (Double.isInfinite(data.getValue(0))) {
                return null;
            }
            myData.E(data);
            myData.TE(1.0/numAtoms);
            return myData;
        }

        protected IEtomicaDataInfo processDataInfo(IEtomicaDataInfo inputDataInfo) {
            return inputDataInfo;
        }
        
        public int numAtoms;
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
            int dim = 3;
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
