package etomica.modules.ljmd;
import java.awt.BorderLayout;

import javax.swing.JPanel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import etomica.action.Action;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorHistory;
import etomica.data.DataFork;
import etomica.data.DataPump;
import etomica.data.DataSink;
import etomica.data.DataSourceCountTime;
import etomica.data.DataSourceFunction;
import etomica.data.DataSourceRmsVelocity;
import etomica.data.DataSourceUniform;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.DataSourceUniform.LimitType;
import etomica.data.meter.MeterDensity;
import etomica.data.meter.MeterEnergy;
import etomica.data.meter.MeterKineticEnergy;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPressure;
import etomica.data.meter.MeterRDF;
import etomica.data.meter.MeterTemperature;
import etomica.data.types.DataFunction;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.ConstantsGraphic;
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
    
    public LjmdGraphic(final Ljmd sim) {
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
            public String getLabel() {return "Phase";}
        }));
        
   /*     DisplayTimer timer = new DisplayTimer(integrator);
        timer.setLabelType(DisplayBox.BORDER);
        timer.setUnit(new Unit(LennardJones.Time.UNIT));
	*/    
	    //meters and displays
        MeterRDF rdfMeter = new MeterRDF(sim.space);
        rdfMeter.getXDataSource().setXMax(4.0);
        rdfMeter.setPhase(sim.phase);
        AccumulatorAverage rdfAverage = new AccumulatorAverage(10);
        DataPump rdfPump = new DataPump(rdfMeter,rdfAverage);
        IntervalActionAdapter rdfAdapter =  new IntervalActionAdapter(rdfPump);
        rdfAdapter.setActionInterval(20);
        sim.integrator.addListener(rdfAdapter);
        sim.register(rdfMeter,rdfPump);
        
        DisplayPlot rdfPlot = new DisplayPlot();
        rdfAverage.addDataSink(rdfPlot.getDataSet(),new StatType[]{StatType.AVERAGE});
        rdfAverage.setPushInterval(1);
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
        rmsAverage.addDataSink(vPlot.getDataSet(), new StatType[]{StatType.AVERAGE});
        vPlot.setDoLegend(false);
//      vPlot.getPlot().setYRange(0.0,0.8);
        vPlot.getPlot().setTitle("Velocity distribution");
        vPlot.setDoLegend(true);
		
		final MaxwellBoltzmann.Distribution mbDistribution = new MaxwellBoltzmann.Distribution(sim.space,sim.integrator.getTemperature(),((AtomTypeLeaf)sim.species.getMoleculeType()).mass);
		final DataSourceFunction mbSource = new DataSourceFunction("Maxwell Boltzmann Distribution",
                Null.DIMENSION, mbDistribution, 100, "Speed", new DimensionRatio(Length.DIMENSION,Time.DIMENSION));
		DataSourceUniform mbX = mbSource.getXSource();
		mbX.setTypeMax(LimitType.HALF_STEP);
		mbX.setTypeMin(LimitType.HALF_STEP);
		mbX.setNValues(((DataFunction.Factory)meterVelocity.getDataInfo().getDataFactory()).getArrayShape()[0]);
		mbX.setXMin(vMin);
		mbX.setXMax(vMax);
		mbSource.update();
        DataPump mbPump = new DataPump(mbSource,vPlot.getDataSet());
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
		
		
/*		double[] xHist = vHistogram.xValues();
		double[] xMB = mbSource.xValues();
		for(int i=0; i<xHist.length; i++) {
		    System.out.println(xHist[i]+"  "+xMB[i]);
		}
*/		
        DataSourceCountTime timeCounter = new DataSourceCountTime();
        sim.integrator.addListener(timeCounter);
		
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
		
	    MeterDensity densityMeter = new MeterDensity(sim.space);
        densityMeter.setPhase(sim.phase);
	    DisplayBox densityBox = new DisplayBox();
        DataPump densityPump = new DataPump(densityMeter, densityBox);
        IntervalActionAdapter densityAdapter = new IntervalActionAdapter(densityPump);
        densityAdapter.setActionInterval(10);
        sim.integrator.addListener(densityAdapter);
        sim.register(densityMeter,densityPump);
	    densityBox.setLabel("Number density");
	    
		MeterEnergy eMeter = new MeterEnergy(sim.potentialMaster);
        eMeter.setPhase(sim.phase);
        AccumulatorHistory energyHistory = new AccumulatorHistory();
        energyHistory.setTimeDataSource(timeCounter);
        DataPump energyPump = new DataPump(eMeter, energyHistory);
        IntervalActionAdapter energyAdapter = new IntervalActionAdapter(energyPump);
        energyAdapter.setActionInterval(60);
        energyHistory.setPushInterval(5);
        sim.integrator.addListener(energyAdapter);
        sim.register(eMeter,energyPump);
		
		MeterPotentialEnergy peMeter = new MeterPotentialEnergy(sim.potentialMaster);
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
        energyHistory.setDataSink(ePlot.getDataSet());
        peHistory.setDataSink(ePlot.getDataSet());
        keHistory.setDataSink(ePlot.getDataSet());
		
		ePlot.setDoLegend(true);
		
		MeterPressure pMeter = new MeterPressure(sim.space);
        pMeter.setIntegrator(sim.integrator);
        pMeter.setIncludeLrc(true);
        AccumulatorAverage pAccumulator = new AccumulatorAverage(sim);
        DataPump pPump = new DataPump(pMeter,pAccumulator);
        IntervalActionAdapter pAdapter = new IntervalActionAdapter(pPump);
        pAdapter.setActionInterval(60);
        pAccumulator.setPushInterval(10);
        sim.register(pMeter,pPump);
        sim.integrator.addListener(pAdapter);
        
        DisplayBoxesCAE pDisplay = new DisplayBoxesCAE();
        pDisplay.setAccumulator(pAccumulator);
        DisplayBoxesCAE peDisplay = new DisplayBoxesCAE();
        peDisplay.setAccumulator(peAccumulator);

        final DeviceNSelector nSlider = new DeviceNSelector(sim,sim.species.getAgent(sim.phase));
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
        bigPanel.add(displayPhase.graphic());
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
        
   //     JPanel timePanel = (JPanel)timer.graphic();
   //     timePanel.setBorder(new javax.swing.border.TitledBorder("Time (LJ)"));
  //      timePanel.add(timer.graphic());
        
 //   	final javax.swing.JTabbedPane sliderPanel = new javax.swing.JTabbedPane();
 //       sliderPanel.add("Number of atoms",nSlider.graphic());
        //panel for all the controls
//        JPanel controlPanel = new JPanel(new java.awt.GridLayout(0,1));
        JPanel controlPanel = new JPanel(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gbc2 = new java.awt.GridBagConstraints();
//        gbc2.gridx = 0;
//        gbc2.gridy = java.awt.GridBagConstraints.RELATIVE;
        gbc2.gridy = 0;
        gbc2.gridx = java.awt.GridBagConstraints.RELATIVE;
        controlPanel.add(startPanel,gbc2);
    //    controlPanel.add(timer.graphic(),gbc2);
    //    controlPanel.add(timePanel,gbc2);
        controlPanel.add(temperaturePanel, gbc2);
        controlPanel.add(nSlider.graphic(), gbc2);
 //       controlPanel.add(sliderPanel, gbc2);
        controlPanel.add(densityBox.graphic(), gbc2);
 //       controlPanel.add(pressureBox.graphic(), gbc2);
        
        panel.add(controlPanel, java.awt.BorderLayout.NORTH);
        panel.add(bigPanel, java.awt.BorderLayout.EAST);
//        panel().add(displayPanel, java.awt.BorderLayout.EAST);
		
		//***************set all the colors******************
		java.awt.Color background = ConstantsGraphic.KHAKI.brighter().brighter();
		java.awt.Color contrast = ConstantsGraphic.DARK_RED;
		panel.setBackground(background);
		bigPanel.setBackground(background);

		controlPanel.setBackground(background);

		startPanel.setBackground(contrast); //border color
		startPanel.setOpaque(false);

		temperaturePanel.setBackground(contrast); //border color
		temperaturePanel.setOpaque(false);
		tBox.graphic(null).setBackground(background);
		tSelect.graphic(null).setBackground(background);
		densityBox.graphic().setBackground(background);
		nSlider.graphic().setBackground(background);
		
    }//end of constructor    
    
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
            
        Ljmd sim = new Ljmd(space);

        LjmdGraphic ljmdGraphic = new LjmdGraphic(sim);
		SimulationGraphic.makeAndDisplayFrame(ljmdGraphic.panel);
    }//end of main
    
    public static class Applet extends javax.swing.JApplet {

	    public void init() {
	        getRootPane().putClientProperty(
	                        "defeatSystemEventQueueCheck", Boolean.TRUE);
		    getContentPane().add(new LjmdGraphic(new Ljmd(Space2D.getInstance())).panel);
	    }
    }//end of Applet
    
    protected JPanel panel;
}//end of LJMD_Module class


