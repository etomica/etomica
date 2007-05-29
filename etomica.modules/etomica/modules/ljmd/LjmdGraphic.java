package etomica.modules.ljmd;
import java.awt.BorderLayout;

import javax.swing.JPanel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import etomica.action.Action;
import etomica.action.SimulationRestart;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.AccumulatorAverage;
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
import etomica.graphics.DefaultToolbar;
import etomica.graphics.DeviceButton;
import etomica.graphics.DeviceNSelector;
import etomica.graphics.DeviceThermoSelector;
import etomica.graphics.DeviceTrioControllerButton;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxesCAE;
import etomica.graphics.DisplayPhase;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntervalActionAdapter;
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

public class LjmdGraphic {
    
    private boolean doConfigButton;
    
    public void setDoConfigButton(boolean newDoConfigButton) {
        doConfigButton = newDoConfigButton;
    }

    public void init(final Ljmd sim) {
        LJ unitSystem = new LJ();
        Unit tUnit = Energy.DIMENSION.getUnit(unitSystem);

        sim.getDefaults().blockSize = 100;
        sim.getDefaults().doSleep = true;
        sim.activityIntegrate.setDoSleep(true);
        
        sim.register(sim.integrator);
        
        panel = new JPanel();
        panel.setLayout(new BorderLayout());
        
        DeviceTrioControllerButton control = new DeviceTrioControllerButton(sim);
        
        //temperature selector
	    DeviceThermoSelector tSelect = new DeviceThermoSelector(sim,sim.integrator);
//	    tSelect.setTemperatures(new double[] {50.,100.,300.,600.,1000.});
	    tSelect.setTemperatures(new double[] {0.3,0.5,0.7,1.0,1.3,2.0,3.0});
	    tSelect.setUnit(tUnit);
	    tSelect.setSelected(0); //sets adiabatic as selected temperature
	    tSelect.getLabel().setText("Set value");
	    
	    //display of phase, timer
	    final DisplayPhase displayPhase = new DisplayPhase(sim.phase,sim.getDefaults().pixelUnit);
        ColorSchemeByType colorScheme = new ColorSchemeByType();
        colorScheme.setColor(sim.species.getMoleculeType(),java.awt.Color.red);
        displayPhase.setColorScheme(new ColorSchemeByType());
        sim.integrator.addListener(new IntervalActionAdapter(new Action() {
            public void actionPerformed() {displayPhase.repaint();}
        }));
        
   /*     DisplayTimer timer = new DisplayTimer(integrator);
        timer.setLabelType(DisplayBox.BORDER);
        timer.setUnit(new Unit(LennardJones.Time.UNIT));
	*/    
	    //meters and displays
        final MeterRDF rdfMeter = new MeterRDF(sim.getSpace());
        IntervalActionAdapter rdfIAA = new IntervalActionAdapter(new Action() {
            public void actionPerformed() {rdfMeter.actionPerformed();}
        });
        rdfIAA.setActionInterval(10);
        sim.integrator.addListener(rdfIAA);
        rdfMeter.getXDataSource().setXMax(4.0);
        rdfMeter.setPhase(sim.phase);
        DisplayPlot rdfPlot = new DisplayPlot();
        DataPump rdfPump = new DataPump(rdfMeter,rdfPlot.getDataSet().makeDataSink());
        IntervalActionAdapter rdfAdapter =  new IntervalActionAdapter(rdfPump);
        rdfAdapter.setActionInterval(10);
        sim.integrator.addListener(rdfAdapter);
        sim.register(rdfMeter,rdfPump);
        
        rdfPlot.setDoLegend(false);
        rdfPlot.getPlot().setTitle("Radial Distribution Function");
		
		//add meter and display for current kinetic temperature
		sim.getDefaults().historyPeriod = 1000;

		
		//velocity distribution
        double vMin = 0;
        double vMax = 4;
		DataSourceRmsVelocity meterVelocity = new DataSourceRmsVelocity(new HistogramSimple(100,new DoubleRange(0,4)));
        meterVelocity.setIterator(new AtomIteratorLeafAtoms(sim.phase));
        AccumulatorAverage rmsAverage = new AccumulatorAverage(10);
        DataPump velocityPump = new DataPump(meterVelocity, rmsAverage);
        IntervalActionAdapter velocityAdapter = new IntervalActionAdapter(velocityPump);
        rmsAverage.setPushInterval(10);
        sim.integrator.addListener(velocityAdapter);
        velocityAdapter.setActionInterval(10);
        sim.register(meterVelocity,velocityPump);
        
        final DisplayPlot vPlot = new DisplayPlot();
        rmsAverage.addDataSink(vPlot.getDataSet().makeDataSink(), new StatType[]{StatType.AVERAGE});
        vPlot.setDoLegend(false);
//      vPlot.getPlot().setYRange(0.0,0.8);
        vPlot.getPlot().setTitle("Velocity distribution");
        vPlot.setDoLegend(true);
		
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
        IntervalActionAdapter mbAdapter = new IntervalActionAdapter(mbPump);
        mbAdapter.setActionInterval(100);
        sim.integrator.addListener(mbAdapter);

		//listner to update Maxwell-Boltzmann plot when temperature is changed
        tSelect.getSelector().addItemListener(new java.awt.event.ItemListener() {
		    public void itemStateChanged(java.awt.event.ItemEvent event) {
		        Object item = event.getItem();
		        if(item instanceof Double) {
		            mbDistribution.setTemperature(((Double)item).doubleValue());
		            mbSource.update();
		            vPlot.doUpdate();
		            vPlot.repaint();
		        }
		    }
		});//end of addItemListener
		
        control.getReinitButton().setPostAction(new Action() {
            public void actionPerformed() {
                displayPhase.repaint();
                rdfMeter.reset();
            }
        });

        control.getResetAveragesButton().setPostAction(new Action() {
            public void actionPerformed() {
                rdfMeter.reset();
            }
        });

        
/*		double[] xHist = vHistogram.xValues();
		double[] xMB = mbSource.xValues();
		for(int i=0; i<xHist.length; i++) {
		    System.out.println(xHist[i]+"  "+xMB[i]);
		}
*/		
        DataSourceCountTime timeCounter = new DataSourceCountTime(sim.integrator);
		
		MeterTemperature thermometer = new MeterTemperature();
        thermometer.setPhase(sim.phase);
        DataFork temperatureFork = new DataFork();
        DataPump temperaturePump = new DataPump(thermometer,temperatureFork);
        IntervalActionAdapter temperatureAdapter = new IntervalActionAdapter(temperaturePump);
        temperatureAdapter.setActionInterval(10);
        sim.integrator.addListener(temperatureAdapter);
        sim.register(meterVelocity,velocityPump);
        AccumulatorHistory temperatureHistory = new AccumulatorHistory();
        temperatureHistory.setTimeDataSource(timeCounter);
		DisplayBox tBox = new DisplayBox();
		temperatureFork.setDataSinks(new DataSink[]{tBox,temperatureHistory});
		tBox.setUnit(tUnit);
		tBox.setLabel("Measured value");
		tBox.setLabelPosition(CompassDirection.NORTH);
		
	    MeterDensity densityMeter = new MeterDensity(sim.getSpace());
        densityMeter.setPhase(sim.phase);
	    DisplayBox densityBox = new DisplayBox();
        DataPump densityPump = new DataPump(densityMeter, densityBox);
        IntervalActionAdapter densityAdapter = new IntervalActionAdapter(densityPump);
        densityAdapter.setActionInterval(10);
        sim.integrator.addListener(densityAdapter);
        sim.register(densityMeter,densityPump);
	    densityBox.setLabel("Number density");
	    
		MeterEnergy eMeter = new MeterEnergy(sim.getPotentialMaster());
        eMeter.setPhase(sim.phase);
        AccumulatorHistory energyHistory = new AccumulatorHistory();
        energyHistory.setTimeDataSource(timeCounter);
        DataPump energyPump = new DataPump(eMeter, energyHistory);
        IntervalActionAdapter energyAdapter = new IntervalActionAdapter(energyPump);
        energyAdapter.setActionInterval(60);
        energyHistory.setPushInterval(5);
        sim.integrator.addListener(energyAdapter);
        sim.register(eMeter,energyPump);
		
		MeterPotentialEnergy peMeter = new MeterPotentialEnergy(sim.getPotentialMaster());
        peMeter.setPhase(sim.phase);
        AccumulatorHistory peHistory = new AccumulatorHistory();
        peHistory.setTimeDataSource(timeCounter);
        AccumulatorAverage peAccumulator = new AccumulatorAverage(sim);
        peAccumulator.setPushInterval(10);
        DataFork peFork = new DataFork(new DataSink[]{peHistory, peAccumulator});
        DataPump pePump = new DataPump(peMeter, peFork);
        IntervalActionAdapter peAdapter = new IntervalActionAdapter(pePump);
        peAdapter.setActionInterval(60);
        peHistory.setPushInterval(5);
        sim.register(peMeter,pePump);
        sim.integrator.addListener(peAdapter);
		
		MeterKineticEnergy keMeter = new MeterKineticEnergy();
        keMeter.setPhase(sim.phase);
        AccumulatorHistory keHistory = new AccumulatorHistory();
        keHistory.setTimeDataSource(timeCounter);
        DataPump kePump = new DataPump(keMeter, keHistory);
        IntervalActionAdapter keAdapter = new IntervalActionAdapter(kePump);
        keAdapter.setActionInterval(60);
        keHistory.setPushInterval(5);
        sim.register(keMeter,kePump);
        sim.integrator.addListener(keAdapter);
        
        DisplayPlot ePlot = new DisplayPlot();
        energyHistory.setDataSink(ePlot.getDataSet().makeDataSink());
        ePlot.setLegend(new DataTag[]{energyHistory.getTag()}, "Total");
        peHistory.setDataSink(ePlot.getDataSet().makeDataSink());
        ePlot.setLegend(new DataTag[]{peHistory.getTag()}, "Potential");
        keHistory.setDataSink(ePlot.getDataSet().makeDataSink());
        ePlot.setLegend(new DataTag[]{keHistory.getTag()}, "Kinetic");
		
		ePlot.setDoLegend(true);
		
        MeterPressureTensorFromIntegrator pMeter = new MeterPressureTensorFromIntegrator();
        pMeter.setIntegrator(sim.integrator);
        AccumulatorAverage pAccumulator = new AccumulatorAverage(sim);
        DataProcessorTensorTrace tracer = new DataProcessorTensorTrace();
        DataPump pPump = new DataPump(pMeter, tracer);
        tracer.setDataSink(pAccumulator);
        IntervalActionAdapter pAdapter = new IntervalActionAdapter(pPump);
        pAccumulator.setPushInterval(10);
        sim.register(pMeter,pPump);
        sim.integrator.addListener(pAdapter);

        DisplayBoxesCAE pDisplay = new DisplayBoxesCAE();
        pDisplay.setAccumulator(pAccumulator);
        DisplayBoxesCAE peDisplay = new DisplayBoxesCAE();
        peDisplay.setAccumulator(peAccumulator);

        final DeviceNSelector nSlider = new DeviceNSelector(sim.getController());
        nSlider.setResetAction(new SimulationRestart(sim));
        nSlider.setSpeciesAgent(sim.species.getAgent(sim.phase));
        nSlider.setMinimum(1);
        nSlider.setMaximum(224);
        nSlider.setLabel("Number of atoms");
        nSlider.setShowBorder(true);
        // add a listener to adjust the thermostat interval for different
        // system sizes (since we're using ANDERSEN_SINGLE.  Smaller systems 
        // don't need as much thermostating.
        nSlider.getSlider().addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent evt) {
                int n = (int)nSlider.getValue();
                sim.integrator.setThermostatInterval(400/n);
                displayPhase.repaint();
            }
        });

        //************* Lay out components ****************//
        
        //tabbed pane for the big displays
        javax.swing.JPanel bigPanel = new javax.swing.JPanel();
        
    	final javax.swing.JTabbedPane displayPanel = new javax.swing.JTabbedPane();
    	try {
        	((javax.swing.JComponent)displayPhase.graphic()).setPreferredSize(new java.awt.Dimension(300,300));
        } catch(ClassCastException e) {}
        displayPhase.setScale(0.7);
        bigPanel.add(displayPanel);
        
//    	javax.swing.JPanel displayPhasePanel = new javax.swing.JPanel();//3d
//    	displayPhasePanel.add(displayPhase.graphic());//3d
//    	displayPanel.add(displayPhase.getLabel(), displayPhasePanel);//3d
        javax.swing.JPanel thermoPanel = new javax.swing.JPanel();
        thermoPanel.add(pDisplay.graphic());
        thermoPanel.add(peDisplay.graphic());
        
        displayPanel.add("Thermo",thermoPanel);
    	displayPanel.add("RDF", rdfPlot.graphic());
    	displayPanel.add("Velocity",vPlot.graphic());
    	
    	javax.swing.JPanel energyPanel = new javax.swing.JPanel(new java.awt.GridLayout(0,1));
 /*   	pePlot.getPlot().setSize(300,200);
    	kePlot.getPlot().setSize(300,200);
    	pePlot.getPlot().setTopPadding(0);
    	kePlot.getPlot().setTopPadding(0);
    	pePlot.getPlot().setTitle("Potential");
    	kePlot.getPlot().setTitle("Kinetic");
 */   //	energyPanel.add(pePlot.graphic());
    //	energyPanel.add(kePlot.graphic());
        ePlot.getPlot().setTitle("Energy History");
        energyPanel.add(ePlot.graphic());
    	displayPanel.add("Energy",energyPanel);
//        displayPanel.add(plot.getLabel(), plot.graphic());
        
        
        displayPanel.addChangeListener(
            new javax.swing.event.ChangeListener() {
                public void stateChanged(javax.swing.event.ChangeEvent event) {
                    displayPanel.invalidate();
                    displayPanel.validate();
                }
        });
        
        //panel for the start buttons
        JPanel startPanel = (JPanel)control.graphic();
        
        //panel for the temperature control/display
        JPanel temperaturePanel = new JPanel(new java.awt.GridBagLayout());
//        temperaturePanel.setBorder(new javax.swing.border.TitledBorder("Temperature (K)"));
        temperaturePanel.setBorder(new javax.swing.border.TitledBorder("Temperature (\u03B5)"));
        java.awt.GridBagConstraints gbc1 = new java.awt.GridBagConstraints();
        gbc1.gridx = 0;  gbc1.gridy = 1;
        gbc1.gridwidth = 1;
        temperaturePanel.add(tSelect.graphic(null),gbc1);
        gbc1.gridx = 1;  gbc1.gridy = 1;
        temperaturePanel.add(tBox.graphic(null),gbc1);
        
        JPanel controlPanel = new JPanel(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gbc2 = new java.awt.GridBagConstraints();
        gbc2.gridx = 0;
        gbc2.gridy = java.awt.GridBagConstraints.RELATIVE;
        controlPanel.add(startPanel,gbc2);
        controlPanel.add(temperaturePanel, gbc2);
        controlPanel.add(nSlider.graphic(), gbc2);
        JPanel anotherPanel = new JPanel(new java.awt.GridBagLayout());
        controlPanel.add(anotherPanel, gbc2);
        gbc2.gridx = 0;
        gbc2.gridy = java.awt.GridBagConstraints.RELATIVE;
        anotherPanel.add(densityBox.graphic(), gbc2);
        
        if (doConfigButton) {
            DeviceButton configButton = new DeviceButton(sim.getController());
            configButton.setLabel("Show Config");
            configButton.setAction(new ActionConfigWindow(sim.phase));
            anotherPanel.add(configButton.graphic(), gbc2);
        }

        DefaultToolbar tb = new DefaultToolbar(panel, "Lennard-Jones Molecular Dynamics");
        panel.add(tb.graphic(), BorderLayout.NORTH);
        JPanel subPanel = new JPanel();
        panel.add(bigPanel, BorderLayout.SOUTH);
        panel.add(controlPanel, BorderLayout.WEST);
        panel.add(displayPhase.graphic());
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
            
        LjmdGraphic ljmdGraphic = new LjmdGraphic();
        ljmdGraphic.setDoConfigButton(true);
        ljmdGraphic.init(new Ljmd(space));
		SimulationGraphic.makeAndDisplayFrame
		        (ljmdGraphic.panel, "Lennard-Jones Molecular Dynamics");
    }
    
    public static class Applet extends javax.swing.JApplet {

        public void init() {
	        getRootPane().putClientProperty(
	                        "defeatSystemEventQueueCheck", Boolean.TRUE);
            LjmdGraphic ljmdGraphic = new LjmdGraphic();
            String doConfigButtonStr = getParameter("doConfigButton");
            if (doConfigButtonStr != null) {
                ljmdGraphic.setDoConfigButton(Boolean.valueOf(doConfigButtonStr).booleanValue());
            }
            ljmdGraphic.init(new Ljmd(Space2D.getInstance()));
		    getContentPane().add(ljmdGraphic.panel);
	    }

        private static final long serialVersionUID = 1L;
    }
    
    protected JPanel panel;
    
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


