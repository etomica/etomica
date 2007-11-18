package etomica.modules.ljmd;

import java.awt.GridBagConstraints;
import java.util.ArrayList;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import etomica.action.Action;
import etomica.action.SimulationRestart;
import etomica.atom.AtomTypeLeaf;
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
import etomica.data.meter.MeterKineticEnergy;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPressureTensorFromIntegrator;
import etomica.data.meter.MeterRDF;
import etomica.data.meter.MeterTemperature;
import etomica.data.types.DataDouble;
import etomica.data.types.DataTensor;
import etomica.data.types.DataFunction.DataInfoFunction;
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

    public LjmdGraphic(final Ljmd simulation) {

    	super(simulation, TABBED_PANE, APP_NAME, REPAINT_INTERVAL);

        ArrayList dataStreamPumps = getController().getDataStreamPumps();
        
    	this.sim = simulation;

        LJ unitSystem = new LJ();
        Unit tUnit = Energy.DIMENSION.getUnit(unitSystem);

        sim.activityIntegrate.setSleepPeriod(1);
       
	    //display of box, timer
        ColorSchemeByType colorScheme = new ColorSchemeByType();
        colorScheme.setColor(sim.species.getMoleculeType(),java.awt.Color.red);
        getDisplayBox(sim.box).setColorScheme(new ColorSchemeByType());
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
		
		final MaxwellBoltzmann.Distribution mbDistribution = new MaxwellBoltzmann.Distribution(sim.getSpace(),sim.integrator.getTemperature(),((AtomTypeLeaf)sim.species.getMoleculeType()).getMass());
		final DataSourceFunction mbSource = new DataSourceFunction("Maxwell Boltzmann Distribution",
                Null.DIMENSION, mbDistribution, 100, "Speed", new DimensionRatio(Length.DIMENSION,Time.DIMENSION));
		DataSourceUniform mbX = mbSource.getXSource();
		mbX.setTypeMax(LimitType.HALF_STEP);
		mbX.setTypeMin(LimitType.HALF_STEP);
		mbX.setNValues(((DataInfoFunction)meterVelocity.getDataInfo()).getLength());
		mbX.setXMin(vMin);
		mbX.setXMax(vMax);
		mbSource.update();
        DataPump mbPump = new DataPump(mbSource,vPlot.getDataSet().makeDataSink());
        sim.integrator.addIntervalAction(mbPump);
        sim.integrator.setActionInterval(mbPump, 100);
		
        getController().getReinitButton().setPostAction(new Action() {
            public void actionPerformed() {
                getDisplayBox(sim.box).repaint();
                rdfMeter.reset();
            }
        });

        getController().getResetAveragesButton().setPostAction(new Action() {
            public void actionPerformed() {
                rdfMeter.reset();
            }
        });

        DataSourceCountTime timeCounter = new DataSourceCountTime(sim.integrator);

        //add meter and display for current kinetic temperature

		MeterTemperature thermometer = new MeterTemperature();
        thermometer.setBox(sim.box);
        DataFork temperatureFork = new DataFork();
        final DataPump temperaturePump = new DataPump(thermometer,temperatureFork);
        sim.integrator.addIntervalAction(temperaturePump);
        sim.integrator.setActionInterval(temperaturePump, 10);
        final AccumulatorHistory temperatureHistory = new AccumulatorHistory();
        temperatureHistory.setTimeDataSource(timeCounter);
		final DisplayTextBox tBox = new DisplayTextBox();
		temperatureFork.setDataSinks(new DataSink[]{tBox,temperatureHistory});
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
	    
		MeterEnergy eMeter = new MeterEnergy(sim.integrator.getPotential());
        eMeter.setBox(sim.box);
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
        DataFork peFork = new DataFork(new DataSink[]{peHistory, peAccumulator});
        DataPump pePump = new DataPump(peMeter, peFork);
        sim.integrator.addIntervalAction(pePump);
        sim.integrator.setActionInterval(pePump, 60);
        peHistory.setPushInterval(5);
        dataStreamPumps.add(pePump);
		
		MeterKineticEnergy keMeter = new MeterKineticEnergy();
        keMeter.setBox(sim.box);
        AccumulatorHistory keHistory = new AccumulatorHistory();
        keHistory.setTimeDataSource(timeCounter);
        DataPump kePump = new DataPump(keMeter, keHistory);
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
		
        MeterPressureTensorFromIntegrator pMeter = new MeterPressureTensorFromIntegrator();
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

        final DeviceNSelector nSlider = new DeviceNSelector(sim.getController());
        nSlider.setResetAction(new SimulationRestart(sim));
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
        nSlider.getSlider().addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent evt) {
                int n = (int)nSlider.getValue();
                if(n == 0) {
                	sim.integrator.setThermostatInterval(400);
                }
                else {
                  sim.integrator.setThermostatInterval(400/n);
                }

                getDisplayBox(sim.box).repaint();
            }
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

	    ChangeListener temperatureListener = new ChangeListener() {
		    public void stateChanged(ChangeEvent event) {

		        mbDistribution.setTemperature(temperatureSelect.getTemperature());
		        mbSource.update();
		        vPlot.doUpdate();
		        vPlot.repaint();
		    }
		};
		temperatureSelect.addTemperatureSliderListener(temperatureListener);

        // show config button
        DeviceButton configButton = new DeviceButton(sim.getController());
        configButton.setLabel("Show Config");
        configButton.setAction(new ActionConfigWindow(sim.box));

        Action resetAction = new Action() {
        	public void actionPerformed() {
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
        Space space = Space2D.getInstance();
        if(args.length != 0) {
            try {
                int D = Integer.parseInt(args[0]);
                if (D == 3) {
                    space = Space3D.getInstance();
                }
            } catch(NumberFormatException e) {}
        }

        LjmdGraphic ljmdGraphic = new LjmdGraphic(new Ljmd(space));
		SimulationGraphic.makeAndDisplayFrame
		        (ljmdGraphic.getPanel(), APP_NAME);
    }
    
    public static class Applet extends javax.swing.JApplet {

        public void init() {
	        getRootPane().putClientProperty(
	                        "defeatSystemEventQueueCheck", Boolean.TRUE);
            LjmdGraphic ljmdGraphic = new LjmdGraphic(new Ljmd(Space2D.getInstance()));

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
        
        protected Data processData(Data inputData) {
            // take the trace and divide by the dimensionality
            data.x = ((DataTensor)inputData).x.trace()/((DataTensor)inputData).x.D();
            return data;
        }

        protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
            dataInfo = new DataDouble.DataInfoDouble(inputDataInfo.getLabel(), inputDataInfo.getDimension());
            return dataInfo;
        }

        public DataPipe getDataCaster(IDataInfo inputDataInfo) {
            if (!(inputDataInfo instanceof DataTensor.DataInfoTensor)) {
                throw new IllegalArgumentException("Gotta be a DataInfoTensor");
            }
            return null;
        }

        private static final long serialVersionUID = 1L;
        protected final DataDouble data;
    }

}


