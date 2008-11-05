package etomica.modules.ljmd;

import java.awt.GridBagConstraints;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;

import etomica.api.IAction;
import etomica.api.IData;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.AccumulatorHistory;
import etomica.data.DataFork;
import etomica.data.DataPipe;
import etomica.data.DataProcessor;
import etomica.data.DataPump;
import etomica.data.DataSourceCountTime;
import etomica.data.DataSourceFunction;
import etomica.data.DataSourceRmsVelocity;
import etomica.data.DataSourceUniform;
import etomica.data.DataTag;
import etomica.data.IDataSink;
import etomica.data.IEtomicaDataInfo;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.DataSourceUniform.LimitType;
import etomica.data.meter.MeterDensity;
import etomica.data.meter.MeterEnergy;
import etomica.data.meter.MeterKineticEnergy;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPressureTensorFromIntegrator;
import etomica.data.meter.MeterRDF;
import etomica.data.meter.MeterTemperature;
import etomica.data.types.DataDouble;
import etomica.data.types.DataTensor;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.exception.ConfigurationOverlapException;
import etomica.graphics.ActionConfigWindow;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DeviceButton;
import etomica.graphics.DeviceNSelector;
import etomica.graphics.DeviceThermoSlider;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.space.Space;
import etomica.space2d.Space2D;
import etomica.space3d.Space3D;
import etomica.statmech.MaxwellBoltzmann;
import etomica.units.DimensionRatio;
import etomica.units.Energy;
import etomica.units.Length;
import etomica.units.Null;
import etomica.units.Time;
import etomica.units.Unit;
import etomica.units.systems.LJ;
import etomica.util.DoubleRange;
import etomica.util.HistogramSimple;
import etomica.util.Constants.CompassDirection;

public class LjmdGraphic extends SimulationGraphic {

    private final static String APP_NAME = "Lennard-Jones Molecular Dynamics";
    private final static int REPAINT_INTERVAL = 20;
    private DeviceThermoSlider temperatureSelect;
    protected Ljmd sim;
    
    private boolean showConfig = false;

    public LjmdGraphic(final Ljmd simulation, Space _space) {

    	super(simulation, TABBED_PANE, APP_NAME, REPAINT_INTERVAL, _space, simulation.getController());

        ArrayList dataStreamPumps = getController().getDataStreamPumps();
        
    	this.sim = simulation;

        LJ unitSystem = new LJ();
        Unit tUnit = Energy.DIMENSION.getUnit(unitSystem);

        sim.activityIntegrate.setSleepPeriod(1);
       
	    //display of box, timer
        ColorSchemeByType colorScheme = new ColorSchemeByType(sim);
        colorScheme.setColor(sim.species.getLeafType(),java.awt.Color.red);
        getDisplayBox(sim.box).setColorScheme(new ColorSchemeByType(sim));
//        sim.integrator.addListener(new IntervalActionAdapter(this.getDisplayBoxPaintAction(sim.box)));

	    //meters and displays
        final MeterRDF rdfMeter = new MeterRDF(sim.getSpace());
        sim.integrator.addIntervalAction(rdfMeter);
        sim.integrator.setActionInterval(rdfMeter, 10);
        rdfMeter.getXDataSource().setXMax(4.0);
        rdfMeter.setBox(sim.box);
        DisplayPlot rdfPlot = new DisplayPlot();
        DataPump rdfPump = new DataPump(rdfMeter,rdfPlot.getDataSet().makeDataSink());
        sim.integrator.addIntervalAction(rdfPump);
        sim.integrator.setActionInterval(rdfPump, 10);
        dataStreamPumps.add(rdfPump);
        
        rdfPlot.setDoLegend(false);
        rdfPlot.getPlot().setTitle("Radial Distribution Function");
        rdfPlot.setLabel("RDF");
        
        final IAction resetDataAction = new IAction(){
            public void actionPerformed() {
                getController().getSimRestart().getDataResetAction();
                rdfMeter.reset();
            }
        };

		//velocity distribution
        double vMin = 0;
        double vMax = 4;
		DataSourceRmsVelocity meterVelocity = new DataSourceRmsVelocity(new HistogramSimple(100,new DoubleRange(0,4)));
        meterVelocity.setIterator(new AtomIteratorLeafAtoms(sim.box));
        AccumulatorAverage rmsAverage = new AccumulatorAverageFixed(10);
        DataPump velocityPump = new DataPump(meterVelocity, rmsAverage);
        sim.integrator.addIntervalAction(velocityPump);
        sim.integrator.setActionInterval(velocityPump, 10);
        rmsAverage.setPushInterval(10);
        dataStreamPumps.add(velocityPump);
        
        final DisplayPlot vPlot = new DisplayPlot();
        rmsAverage.addDataSink(vPlot.getDataSet().makeDataSink(), new StatType[]{StatType.AVERAGE});
        vPlot.setDoLegend(false);
        vPlot.getPlot().setTitle("Velocity Distribution");
        vPlot.setDoLegend(true);
        vPlot.setLabel("Velocity");
		
		final MaxwellBoltzmann.Distribution mbDistribution = new MaxwellBoltzmann.Distribution(sim.getSpace(),sim.integrator.getTemperature(),sim.species.getLeafType().getMass());
		final DataSourceFunction mbSource = new DataSourceFunction("Maxwell Boltzmann Distribution",
                Null.DIMENSION, mbDistribution, 100, "Speed", new DimensionRatio(Length.DIMENSION,Time.DIMENSION));
		DataSourceUniform mbX = mbSource.getXSource();
		mbX.setTypeMax(LimitType.HALF_STEP);
		mbX.setTypeMin(LimitType.HALF_STEP);
		mbX.setNValues(((DataInfoFunction)meterVelocity.getDataInfo()).getLength());
		mbX.setXMin(vMin);
		mbX.setXMax(vMax);
		mbSource.update();
        final DataPump mbPump = new DataPump(mbSource,vPlot.getDataSet().makeDataSink());
		
        DataSourceCountTime timeCounter = new DataSourceCountTime(sim.integrator);

        //add meter and display for current kinetic temperature

		MeterTemperature thermometer = new MeterTemperature(sim.box, space.D());
        DataFork temperatureFork = new DataFork();
        final DataPump temperaturePump = new DataPump(thermometer,temperatureFork);
        sim.integrator.addIntervalAction(temperaturePump);
        sim.integrator.setActionInterval(temperaturePump, 10);
        final AccumulatorHistory temperatureHistory = new AccumulatorHistory();
        temperatureHistory.setTimeDataSource(timeCounter);
		final DisplayTextBox tBox = new DisplayTextBox();
		temperatureFork.setDataSinks(new IDataSink[]{tBox,temperatureHistory});
        tBox.setUnit(tUnit);
		tBox.setLabel("Measured Temperature");
		tBox.setLabelPosition(CompassDirection.NORTH);

		dataStreamPumps.add(temperaturePump);
        tBox.setUnit(tUnit);
		tBox.setLabel("Measured Temperature");
		tBox.setLabelPosition(CompassDirection.NORTH);

		// Number density box
	    MeterDensity densityMeter = new MeterDensity(sim.getSpace());
        densityMeter.setBox(sim.box);
	    final DisplayTextBox densityBox = new DisplayTextBox();
        final DataPump densityPump = new DataPump(densityMeter, densityBox);
        sim.integrator.addIntervalAction(densityPump);
        sim.integrator.setActionInterval(densityPump, 10);
        dataStreamPumps.add(densityPump);
	    densityBox.setLabel("Number Density");
	    
		MeterEnergy eMeter = new MeterEnergy(sim.integrator.getPotential(), sim.box);
        AccumulatorHistory energyHistory = new AccumulatorHistory();
        energyHistory.setTimeDataSource(timeCounter);
        DataPump energyPump = new DataPump(eMeter, energyHistory);
        sim.integrator.addIntervalAction(energyPump);
        sim.integrator.setActionInterval(energyPump, 60);
        energyHistory.setPushInterval(5);
        dataStreamPumps.add(energyPump);
		
		MeterPotentialEnergy peMeter = new MeterPotentialEnergy(sim.integrator.getPotential());
        peMeter.setBox(sim.box);
        AccumulatorHistory peHistory = new AccumulatorHistory();
        peHistory.setTimeDataSource(timeCounter);
        final AccumulatorAverageCollapsing peAccumulator = new AccumulatorAverageCollapsing();
        peAccumulator.setPushInterval(10);
        DataFork peFork = new DataFork(new IDataSink[]{peHistory, peAccumulator});
        DataPump pePump = new DataPump(peMeter, peFork);
        sim.integrator.addIntervalAction(pePump);
        sim.integrator.setActionInterval(pePump, 60);
        peHistory.setPushInterval(5);
        dataStreamPumps.add(pePump);
		
		MeterKineticEnergy keMeter = new MeterKineticEnergy();
        keMeter.setBox(sim.box);
        AccumulatorHistory keHistory = new AccumulatorHistory();
        keHistory.setTimeDataSource(timeCounter);
        DataFork keFork = new DataFork();
        DataPump kePump = new DataPump(keMeter, keFork);
        keFork.addDataSink(keHistory);
        final AccumulatorAverage keAvg = new AccumulatorAverageCollapsing();
        keFork.addDataSink(keAvg);
        sim.integrator.addIntervalAction(kePump);
        sim.integrator.setActionInterval(kePump, 60);
        keHistory.setPushInterval(5);
        dataStreamPumps.add(kePump);
        
        DisplayPlot ePlot = new DisplayPlot();
        energyHistory.setDataSink(ePlot.getDataSet().makeDataSink());
        ePlot.setLegend(new DataTag[]{energyHistory.getTag()}, "Total");
        peHistory.setDataSink(ePlot.getDataSet().makeDataSink());
        ePlot.setLegend(new DataTag[]{peHistory.getTag()}, "Potential");
        keHistory.setDataSink(ePlot.getDataSet().makeDataSink());
        ePlot.setLegend(new DataTag[]{keHistory.getTag()}, "Kinetic");

        ePlot.getPlot().setTitle("Energy History");
		ePlot.setDoLegend(true);
		ePlot.setLabel("Energy");
		
        MeterPressureTensorFromIntegrator pMeter = new MeterPressureTensorFromIntegrator(space);
        pMeter.setIntegrator(sim.integrator);
        final AccumulatorAverageCollapsing pAccumulator = new AccumulatorAverageCollapsing();
        DataProcessorTensorTrace tracer = new DataProcessorTensorTrace();
        final DataPump pPump = new DataPump(pMeter, tracer);
        tracer.setDataSink(pAccumulator);
        sim.integrator.addIntervalAction(pPump);
        pAccumulator.setPushInterval(10);
        dataStreamPumps.add(pPump);

        final DisplayTextBoxesCAE pDisplay = new DisplayTextBoxesCAE();
        pDisplay.setAccumulator(pAccumulator);
        final DisplayTextBoxesCAE peDisplay = new DisplayTextBoxesCAE();
        peDisplay.setAccumulator(peAccumulator);
        
        IAction mbAction = new IAction() {
            public void actionPerformed() {
                if (!sim.integrator.isIsothermal()) {
                    // in adiabatic mode, we want the average kinetic temperature, which we'll
                    // steal from our kinetic energy accumulator
                    double temp = keAvg.getData().getValue(AccumulatorAverage.StatType.AVERAGE.index);
                    temp *= 2.0 / (sim.box.getLeafList().getAtomCount() * space.D());
                    mbDistribution.setTemperature(temp);
                    mbSource.update();
                }
                mbPump.actionPerformed();
            }
        };
            
        sim.integrator.addIntervalAction(mbAction);
        sim.integrator.setActionInterval(mbAction, 100);

        final DeviceNSelector nSlider = new DeviceNSelector(sim.getController());
        nSlider.setSpecies(sim.species);
        nSlider.setBox(sim.box);
        nSlider.setMinimum(0);
        nSlider.setMaximum(225);
        nSlider.setLabel("Number of Atoms");
        nSlider.setShowBorder(true);
        nSlider.setShowValues(true);
        // add a listener to adjust the thermostat interval for different
        // system sizes (since we're using ANDERSEN_SINGLE.  Smaller systems 
        // don't need as much thermostating.
        nSlider.setPostAction(new IAction() {
            public void actionPerformed() {
                int n = (int)nSlider.getValue();
                if(n == 0) {
                	sim.integrator.setThermostatInterval(400);
                }
                else {
                  sim.integrator.setThermostatInterval(400/n);
                }
                
                if (oldN < n) {
                    config.initializeCoordinates(sim.box);
                }
                oldN = n;
                try {
                    sim.integrator.reset();
                }
                catch (ConfigurationOverlapException e) {
                    throw new RuntimeException(e);
                }
                resetDataAction.actionPerformed();
                getDisplayBox(sim.box).repaint();
            }
            
            ConfigurationLattice config = new ConfigurationLattice((space.D() == 2) ? new LatticeOrthorhombicHexagonal(space) : new LatticeCubicFcc(space), space);
            int oldN = sim.box.getMoleculeList().getAtomCount();
        });

        //************* Lay out components ****************//

        GridBagConstraints vertGBC = SimulationPanel.getVertGBC();

        getDisplayBox(sim.box).setScale(0.7);

        //temperature selector
        temperatureSelect = new DeviceThermoSlider(sim.getController());
        temperatureSelect.setPrecision(1);
        temperatureSelect.setMinimum(0.0);
        temperatureSelect.setMaximum(3.0);
        temperatureSelect.setSliderMajorValues(3);
	    temperatureSelect.setUnit(tUnit);
	    temperatureSelect.setAdiabatic();
	    temperatureSelect.setIntegrator(sim.integrator);

        final IAction temperatureAction = new IAction() {
            public void actionPerformed() {
                resetDataAction.actionPerformed();
		        mbDistribution.setTemperature(temperatureSelect.getTemperature());
		        mbSource.update();
		        vPlot.doUpdate();
		        vPlot.repaint();
		    }
		};
		ActionListener isothermalListener = new ActionListener() {
		    public void actionPerformed(ActionEvent event) {
		        // we can't tell if we're isothermal here...  :(
		        // if we're adiabatic, we'll re-set the temperature elsewhere
		        temperatureAction.actionPerformed();
		    }
		};

		temperatureSelect.setSliderPostAction(temperatureAction);
        temperatureSelect.addRadioGroupActionListener(isothermalListener);

        // show config button
        DeviceButton configButton = new DeviceButton(sim.getController());
        configButton.setLabel("Show Config");
        configButton.setAction(new ActionConfigWindow(sim.box));

        IAction resetAction = new IAction() {
        	public void actionPerformed() {
                rdfMeter.reset();

                // Reset density (Density is set and won't change, but
        		// do this anyway)
        		densityPump.actionPerformed();
        		densityBox.repaint();

        		// Reset temperature (THIS IS NOT WORKING)
                temperaturePump.actionPerformed();
//                tBox.putData(temperatureHistory.getData());
                tBox.repaint();

                // IS THIS WORKING?
                pPump.actionPerformed();
                pDisplay.putData(pAccumulator.getData());
                pDisplay.repaint();
                peDisplay.putData(peAccumulator.getData());
                peDisplay.repaint();

        		getDisplayBox(sim.box).graphic().repaint();
        	}
        };

        this.getController().getReinitButton().setPostAction(resetAction);
        this.getController().getResetAveragesButton().setPostAction(resetAction);

        getPanel().controlPanel.add(temperatureSelect.graphic(), vertGBC);
        add(nSlider);
        if(showConfig == true) {
            add(configButton);
        }

    	add(rdfPlot);
    	add(vPlot);
    	add(ePlot);
    	add(densityBox);
    	add(tBox);
    	add(pDisplay);
    	add(peDisplay);

    }

    public static void main(String[] args) {
        Space sp = null;
        if(args.length != 0) {
            try {
                int D = Integer.parseInt(args[0]);
                if (D == 3) {
                    sp = Space3D.getInstance();
                }
                else {
                	sp = Space2D.getInstance();
                }
            } catch(NumberFormatException e) {}
        }
        else {
        	sp = Space2D.getInstance();
        }

        LjmdGraphic ljmdGraphic = new LjmdGraphic(new Ljmd(sp), sp);
		SimulationGraphic.makeAndDisplayFrame
		        (ljmdGraphic.getPanel(), APP_NAME);
    }
    
    public static class Applet extends javax.swing.JApplet {

        public void init() {
	        getRootPane().putClientProperty(
	                        "defeatSystemEventQueueCheck", Boolean.TRUE);
	        Space sp = Space2D.getInstance();
            LjmdGraphic ljmdGraphic = new LjmdGraphic(new Ljmd(sp), sp);

		    getContentPane().add(ljmdGraphic.getPanel());
	    }

        private static final long serialVersionUID = 1L;
    }
    
    /**
     * Inner class to find the total pressure of the system from the pressure
     * tensor.
     */
    public static class DataProcessorTensorTrace extends DataProcessor {

        public DataProcessorTensorTrace() {
            data = new DataDouble();
        }
        
        protected IData processData(IData inputData) {
            // take the trace and divide by the dimensionality
            data.x = ((DataTensor)inputData).x.trace()/((DataTensor)inputData).x.D();
            return data;
        }

        protected IEtomicaDataInfo processDataInfo(IEtomicaDataInfo inputDataInfo) {
            dataInfo = new DataDouble.DataInfoDouble(inputDataInfo.getLabel(), inputDataInfo.getDimension());
            return dataInfo;
        }

        public DataPipe getDataCaster(IEtomicaDataInfo inputDataInfo) {
            if (!(inputDataInfo instanceof DataTensor.DataInfoTensor)) {
                throw new IllegalArgumentException("Gotta be a DataInfoTensor");
            }
            return null;
        }

        private static final long serialVersionUID = 1L;
        protected final DataDouble data;
    }

}


