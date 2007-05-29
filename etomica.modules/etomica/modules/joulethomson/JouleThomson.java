package etomica.modules.joulethomson;

import java.awt.Color;
import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;

import etomica.action.Action;
import etomica.action.AtomAction;
import etomica.action.AtomActionAdapter;
import etomica.action.ResetAccumulators;
import etomica.action.SimulationDataAction;
import etomica.action.SimulationRestart;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.AtomTypeSphere;
import etomica.atom.IAtom;
import etomica.chem.elements.ElementSimple;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorHistory;
import etomica.data.DataPump;
import etomica.data.DataSource;
import etomica.data.DataSourceScalar;
import etomica.data.DataTag;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.meter.MeterDensity;
import etomica.graphics.AboutBoxWindow;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.ColorSchemeTemperature;
import etomica.graphics.DeviceBox;
import etomica.graphics.DeviceButton;
import etomica.graphics.DeviceControllerButton;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DeviceToggleRadioButtons;
import etomica.graphics.DisplayPhase;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTable;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorGear4NPH;
import etomica.integrator.IntervalActionAdapter;
import etomica.integrator.IntegratorGear4NPH.MeterTPH;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.lattice.SpaceLattice;
import etomica.modifier.ModifierFunctionWrapper;
import etomica.modifier.ModifierGeneral;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.space2d.Space2D;
import etomica.units.Bar;
import etomica.units.CompoundUnit;
import etomica.units.Joule;
import etomica.units.Kelvin;
import etomica.units.Liter;
import etomica.units.Meter;
import etomica.units.Mole;
import etomica.units.Prefix;
import etomica.units.PrefixedUnit;
import etomica.units.Pressure;
import etomica.units.Temperature;
import etomica.units.Unit;
import etomica.units.UnitRatio;
import etomica.units.systems.MKS;

/**
 * Joule Thomson Experiment
 * Ed Maginn [University of Notre Dame], David Kofke [University at Buffalo]
 * Permits study of a system undergoing isothermal or isenthalpic expansion and
 * compression. Isenthalpic expansion is known as  Joule-Thomson process, and 
 * is the basis of many common refrigeration processes.
 */
public class JouleThomson extends SimulationGraphic {
    
    DeviceButton restart;
    DisplayPhase displayPhase1;
    DeviceBox timeBox;
    JouleThomsonSim sim;
    
    public JouleThomson() {this(Space2D.getInstance());}
    public JouleThomson(Space space) {
        super(new JouleThomsonSim(space));
        sim = (JouleThomsonSim)getSimulation();
        
        sim.register(sim.integratorJT);
        
        final Unit pUnit, dUnit, dadUnit;
        if (space.D() == 2) {
            Unit areaUnit = new MKS().area();
            dUnit = new UnitRatio(Mole.UNIT, 
                                    areaUnit);
            dadUnit = new UnitRatio(new PrefixedUnit(Prefix.DECI, Mole.UNIT),
                                    areaUnit);
            Unit[] units = new Unit[] {Bar.UNIT, new PrefixedUnit(Prefix.NANO, Meter.UNIT)};
            double[] exponents = new double[] {1.0, 1.0};
            pUnit = new CompoundUnit(units, exponents);
        }
        else {
            dUnit = new UnitRatio(Mole.UNIT, Liter.UNIT);
            dadUnit = new UnitRatio(new PrefixedUnit(Prefix.DECI, Mole.UNIT),
                    Liter.UNIT);
            pUnit = Bar.UNIT;
        }

        final Unit tUnit = Kelvin.UNIT;
        Unit hUnit = new UnitRatio(Joule.UNIT, Mole.UNIT);
 
	    displayPhase1 = new DisplayPhase(sim.phase);
        sim.integratorJT.addListener(new IntervalActionAdapter(new Action() {
            public void actionPerformed() {
                displayPhase1.repaint();
            }
        }));
	    
	    //colorscheme to color atoms blue to red according to their velocity
	    DeviceSlider scaleSlider = null;
        if(space.D() == 2) 
            displayPhase1.setColorScheme(new ColorSchemeTemperature(100, 500));
        else {
            ColorSchemeByType colorScheme = new ColorSchemeByType();
            colorScheme.setColor(sim.species.getMoleculeType(), Color.blue);
            displayPhase1.setColorScheme(colorScheme);
        }
        ModifierFunctionWrapper scaleModifier = new ModifierFunctionWrapper(displayPhase1, "scale");
        scaleModifier.setFunction(new etomica.util.Function.Linear(0.01, 0.0));
        scaleSlider = new DeviceSlider(sim.getController());
        scaleSlider.setModifier(scaleModifier);
        scaleSlider.getSlider().addChangeListener(new javax.swing.event.ChangeListener() {
            public void stateChanged(javax.swing.event.ChangeEvent evt) {
                displayPhase1.graphic().repaint();
            }
        });
        scaleSlider.setMinimum(10);
        scaleSlider.setMaximum(100);
//	    scaleSlider.getSlider().setSnapToTicks(true);
        scaleSlider.getSlider().setValue(100);
        scaleSlider.getSlider().setMajorTickSpacing(10);
        scaleSlider.getSlider().setMinorTickSpacing(5);
        scaleSlider.getSlider().setOrientation(1);
        scaleSlider.getSlider().setLabelTable(scaleSlider.getSlider().createStandardLabels(10));
        scaleSlider.setLabel("Scale (%)");
        scaleSlider.setShowBorder(true);
	    //slider for scale of display
	    
//		panel().setBackground(java.awt.Color.blue);
		

//		DeviceTrioControllerButton control = new DeviceTrioControllerButton(sim);
  
        DeviceToggleRadioButtons isothermalToggleButton = new DeviceToggleRadioButtons(new IntegratorGear4NPH.EnsembleToggler(sim.integrator), "Ensemble", "NPT", "NPH");

/*        DeviceToggleButton sleepToggleButton = new DeviceToggleButton(
            this, new ModulatorBoolean() {
                    public void setBoolean(boolean b) {integrator.setDoSleep(b);}
                    public boolean getBoolean() {return integrator.isDoSleep();}
                    },
                    "Slow", "Fast");
*/		
		//pressure and temperature sliders
//		if(space.D() == 3) integrator.setDoSleep(false);
		
		final DeviceSlider pSlider = new DeviceSlider(sim.getController());
		pSlider.setUnit(pUnit);
        pSlider.setModifier(new ModifierGeneral(sim.integrator,"targetP"));
		pSlider.setLabel("Setpoint Pressure ("+pUnit.symbol()+")");
		pSlider.setMinimum(0);
		pSlider.setMaximum(200);
		pSlider.setNMajor(5);
		pSlider.setShowBorder(true);
		
		final DeviceSlider tSlider = new DeviceSlider(sim.getController());
        tSlider.setModifier(new ModifierGeneral(sim.integrator,"targetT"));
		tSlider.setUnit(tUnit);
		tSlider.setLabel("Setpoint Temperature ("+tUnit.symbol()+")");
		tSlider.setNMajor(4);
		tSlider.setShowBorder(true);
		
        //restart action and button
        etomica.action.SimulationRestart restartAction = new etomica.action.SimulationRestart(sim);
///        restart = new DeviceButton(restartAction);
        restart = new DeviceButton(sim.getController());
        restart.setAction(restartAction);
        restart.setLabel("Restart");
        
        //reset-averages action and button
        SimulationDataAction resetAveragesAction = new SimulationDataAction(sim, new ResetAccumulators());
        DeviceButton resetAverages = new DeviceButton(sim.getController());
        resetAverages.setAction(resetAveragesAction);
        resetAverages.setLabel("Reset averages");

        //temperature, pressure, enthalpy meters
        //XXX we have no way of using this Meter's output  :(
        MeterTPH properties = new MeterTPH(sim.integrator);
//        etomica.MeterScalar tMeter = properties.allMeters()[0];
//        etomica.MeterScalar pMeter = properties.allMeters()[1];
//        etomica.MeterScalar hMeter = properties.allMeters()[2];
//        properties.setHistorying(true);
//        tMeter.setHistorying(true);
//        tMeter.setPhase(phase);
//        pMeter.setHistorying(true);
//        hMeter.setHistorying(true);
        
        //meter and display for density
        final MeterDensity densityMeter = new MeterDensity(space);
        densityMeter.setPhase(sim.phase);
        AccumulatorHistory densityHistoryAcc = new AccumulatorHistory();
        DataPump densityPump = new DataPump(densityMeter, densityHistoryAcc);
        IntervalActionAdapter adapter = new IntervalActionAdapter(densityPump);
        sim.integratorJT.addListener(adapter);
        sim.integratorJT.setEventInterval(20);
        
        //plot of temperature and density histories
//		History tMeterHistory = tMeter.getHistory();
//		History pMeterHistory = pMeter.getHistory();
//		History hMeterHistory = hMeter.getHistory();
//		densityHistory.setLabel("Density ("+dadUnit.symbol()+")");
//		tMeterHistory.setLabel("Temperature ("+tUnit.symbol()+")");
//		pMeterHistory.setLabel("Pressure ("+pUnit.symbol()+")");
//		hMeterHistory.setLabel("Enthalpy ("+hUnit.symbol()+")");


        sim.activityIntegrate.setDoSleep(true);
        sim.activityIntegrate.setSleepPeriod(1);
        System.out.println("Time step: "+sim.integrator.getTimeStep());

        //set-pressure history
        DataSource targetPressureDataSource = new DataSourceScalar("Set-Pressure",Pressure.DIMENSION) {
            public double getDataAsScalar() {
                return pUnit.toSim(pSlider.getValue());
            }
        };
        AccumulatorHistory targetPressureHistory = new AccumulatorHistory();
        DataPump pump = new DataPump(targetPressureDataSource, targetPressureHistory);
        adapter = new IntervalActionAdapter(pump,sim.integratorJT);

        //set-pressure history
        DataSource targetTemperatureDataSource = new DataSourceScalar("Set-Temperature",Temperature.DIMENSION) {
            public double getDataAsScalar() {
                return tUnit.toSim(tSlider.getValue());
            }
        };
        tSlider.setValue(tSlider.getUnit().fromSim(sim.integrator.getTemperature()));
        AccumulatorHistory targetTemperatureHistory = new AccumulatorHistory();
        pump = new DataPump(targetTemperatureDataSource, targetTemperatureHistory);
        adapter = new IntervalActionAdapter(pump,sim.integratorJT);

        //plot of pressure-density-temperature setpoints and averages
        DisplayPlot plot = new DisplayPlot();
        densityHistoryAcc.addDataSink(plot.getDataSet().makeDataSink());
        //add actual temperature history
        targetTemperatureHistory.addDataSink(plot.getDataSet().makeDataSink());
        // add actual pressure history
        targetPressureHistory.addDataSink(plot.getDataSet().makeDataSink());
        plot.setLabel("PVT");
        plot.getPlot().setYLabel("");
        plot.getPlot().setYRange(0.0, 500.);
        plot.setUnit(new DataTag[]{densityHistoryAcc.getTag()}, dadUnit);
        plot.setUnit(new DataTag[]{targetTemperatureHistory.getTag()}, tUnit);
        plot.setUnit(new DataTag[]{targetPressureHistory.getTag()}, pUnit);

        //plot of enthalpy and PT set points
        DisplayPlot plotH = new DisplayPlot();
        targetTemperatureHistory.addDataSink(plotH.getDataSet().makeDataSink());
        targetPressureHistory.addDataSink(plotH.getDataSet().makeDataSink());
        // add enthalpy history
        plotH.setLabel("Enthalpy");
        plotH.getPlot().setYLabel("");
        plot.setUnit(new DataTag[]{targetTemperatureHistory.getTag()}, tUnit);
        plot.setUnit(new DataTag[]{targetPressureHistory.getTag()}, pUnit);
        
        SpeciesChooser speciesChooser = new SpeciesChooser(this);
	    speciesChooser.setSpecies("Methane");
	    
        DisplayTable displayTable = new DisplayTable();
        AccumulatorAverage densityAverage = new AccumulatorAverage(sim);
        pump = new DataPump(densityMeter, densityAverage);
        adapter = new IntervalActionAdapter(pump);
        sim.integratorJT.addListener(adapter);
        densityAverage.addDataSink(displayTable.getDataTable().makeDataSink(), new StatType[]{StatType.MOST_RECENT, StatType.AVERAGE});
        // add temp meter, pressure meter, enthalpy meter to display table
        displayTable.setUnit(new DataTag[]{densityAverage.getTag()}, dUnit);
        // set units for temp meter, pressure meter, enthalpy meter to display table
        
                //*************************************************//
                //************* Lay out components ****************//
                //*************************************************//
        
        //tabbed pane for the big displays
    	final javax.swing.JTabbedPane displayPanel = new javax.swing.JTabbedPane();
	    javax.swing.JPanel displayPhasePanel = new javax.swing.JPanel(new BorderLayout());
	    displayPhasePanel.add(displayPhase1.graphic(),BorderLayout.CENTER);
	    if(space.D() < 3) displayPhasePanel.add(scaleSlider.getSlider(),BorderLayout.EAST);
//   	displayPanel.add(displayPhase1.getLabel(), displayPhase1.graphic());
    	displayPanel.add(displayPhase1.getLabel(), displayPhasePanel);
    	displayPanel.add(displayTable.getLabel(), displayTable.graphic());
        displayPanel.add(plot.getLabel(), plot.graphic());
        displayPanel.add(plotH.getLabel(), plotH.graphic());
        
        displayPanel.addChangeListener(
            new javax.swing.event.ChangeListener() {
                public void stateChanged(javax.swing.event.ChangeEvent event) {
                    displayPanel.invalidate();
                    displayPanel.validate();
                }
        });
        
        timeBox = new DeviceBox();
        timeBox.setModifier(new ModifierGeneral(sim.integrator, "timeStep"));
        
        //panel for the start buttons
//        JPanel startPanel = new JPanel(new java.awt.FlowLayout());
        DeviceControllerButton startButton = new DeviceControllerButton(sim.getController());
///        etomica.Controller.Button startButton = controller1.getButton();
        JPanel startPanel = new JPanel(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gbc0 = new java.awt.GridBagConstraints();
        startPanel.setBorder(new javax.swing.border.TitledBorder("Control"));
        gbc0.gridx = 0; gbc0.gridy = 0;
 //       startPanel.add(controller1.getButton().graphic(null), gbc0);
        startPanel.add(startButton.graphic(null), gbc0);
        gbc0.gridx = 1; gbc0.gridy = 0;
        startPanel.add(restart.graphic(null), gbc0);
        gbc0.gridx = 0; gbc0.gridy = 1; gbc0.gridwidth = 2;
        startPanel.add(resetAverages.graphic(null), gbc0);
        gbc0.gridx = 0; gbc0.gridy = 2; gbc0.gridwidth = 2;
        startPanel.add(isothermalToggleButton.graphic(null), gbc0);
        gbc0.gridx = 0; gbc0.gridy = 3; gbc0.gridwidth = 2;
        startPanel.add(timeBox.graphic(), gbc0);
//        startPanel.add(sleepToggleButton.graphic(null), gbc0);
        
        //panel for sliders
        JPanel sliderPanel = new JPanel(new java.awt.GridLayout(0,1));
        sliderPanel.add(pSlider.graphic(null));
        sliderPanel.add(tSlider.graphic(null));

        //panel for all the controls
//        JPanel controlPanel = new JPanel(new java.awt.GridLayout(0,1));
        JPanel controlPanel = new JPanel(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gbc2 = new java.awt.GridBagConstraints();
        gbc2.gridx = 0;
        gbc2.gridy = java.awt.GridBagConstraints.RELATIVE;
        
        //species selector
        JPanel potentialPanel = new JPanel();
        potentialPanel.add(speciesChooser);
	    potentialPanel.setBorder(new javax.swing.border.TitledBorder("Species selection"));

        controlPanel.add(startPanel,gbc2);
        controlPanel.add(potentialPanel, gbc2);
        controlPanel.add(sliderPanel, gbc2);
        
        panel().removeAll();
        panel().setLayout(new BorderLayout());
        panel().add(controlPanel, BorderLayout.WEST);
        panel().add(displayPanel, BorderLayout.EAST);
        panel().add(addMenu(panel()), BorderLayout.NORTH);

    }

    private JMenuBar addMenu(final JPanel parent) {
    	JMenuBar mBar = new JMenuBar();
    	JMenu fileMenu = new JMenu("File");
    	JMenuItem exitBtn = new JMenuItem("Exit");
    	exitBtn.addActionListener(new ActionListener() {
    		public void actionPerformed(ActionEvent ev) {
    			System.exit(0);
    		}
    	});
    	fileMenu.add(exitBtn);
    	JMenu helpMenu = new JMenu("Help");
    	JMenuItem aboutBtn = new JMenuItem("About Joule Thomson Experiment");
    	aboutBtn.addActionListener(new ActionListener() {
    		public void actionPerformed(ActionEvent ev) {
    			AboutBoxWindow about =
    				new AboutBoxWindow(parent,
    					               "About Joule Thomson Experiment",
    					               new String[] {"Dr. David A. Kofke", "Dr. Andrew Schultz" },
    					               new String[] {"Colin Tedlock", "Robert Rassler" });
    			about.setVisible(true);
    		}
    	});
    	aboutBtn.setEnabled(true);
    	helpMenu.add(aboutBtn);

    	mBar.add(fileMenu);
    	mBar.add(helpMenu);

    	return(mBar);

    }

    //inner class the defines a drop-down menu to select LJ parameters to mimic
    //several real substances
    //the field "names" must be static to pass to super, and that requires this inner class to be static
    public static class SpeciesChooser extends javax.swing.JComboBox implements java.awt.event.ItemListener {
        
        public static final String[] names = new String[] 
            {"Methane", "Argon", "Krypton", "Nitrogen", "Ideal gas"};
        double[] sigma   = new double[] {
                                         4.010, 
                                         3.499, 
                                         3.846, 
                                         3.694, 
                                         2.};
        double[] epsilon = new double[] {
                                         Kelvin.UNIT.toSim(142.87), 
                                         Kelvin.UNIT.toSim(118.13),
                                         Kelvin.UNIT.toSim(162.74), 
                                         Kelvin.UNIT.toSim(96.26),
                                         Kelvin.UNIT.toSim(100.)};
        double[] mass    = new double[] {
                                         16.0, 
                                         39.95,  
                                         83.80, 
                                         28.0, 
                                         1.0};
        JouleThomsonSim sim;
        JouleThomson simGraphic;
        double currentMass = mass[0];
        double currentEps = epsilon[0];
        double currentSig = sigma[0];
        AtomAction updateMass = new AtomActionAdapter() {
            public void actionPerformed(IAtom a) {((ElementSimple)((AtomTypeLeaf)a.getType()).getElement()).setMass(currentMass);}
        };
        SimulationRestart simRestart;
        
        SpeciesChooser(JouleThomson simGraphic) {
            super(names);
            this.simGraphic = simGraphic;
            setPreferredSize(new java.awt.Dimension(150,30));
            setOpaque(false);
            addItemListener(this);
            sim = (JouleThomsonSim)simGraphic.getSimulation();
            currentEps = sim.potential.getEpsilon();
            currentSig = sim.potential.getSigma();
            sigma[0] = currentSig;
            epsilon[0] = currentEps;
            mass[0] = currentMass;
            simRestart = new SimulationRestart(sim);
        }
        
        public void itemStateChanged(java.awt.event.ItemEvent evt) {
            if(evt.getStateChange() == java.awt.event.ItemEvent.DESELECTED) return;
            String speciesName = (String)evt.getItemSelectable().getSelectedObjects()[0];
            setSpecies(speciesName);
//			simGraphic.timeBox.doUpdate();
            simGraphic.displayPhase1.repaint();
        }
        
        public void setSpecies(String speciesName) {
            for(int i=0; i<names.length; i++) {
                if(speciesName.equals(names[i])) {
                    double tStep = sim.integrator.getTimeStep();
                    System.out.println("old "+tStep);
//                    tStep *= sigma[i]/currentSig*Math.sqrt(mass[i]*currentEps/(currentMass*epsilon[i]));
                    tStep *= sigma[i]/currentSig*Math.sqrt(mass[i]/(currentMass));
                    System.out.println("new "+tStep);
                    sim.integrator.setTimeStep(tStep);
                    sim.integratorNVE.setTimeStep(tStep);
                    currentMass = mass[i];
                    ((ElementSimple)((AtomTypeLeaf)sim.species.getMoleculeType()).getElement()).setMass(currentMass);
                    currentEps = epsilon[i];
                    currentSig = sigma[i];
                    ((AtomTypeSphere)sim.species.getMoleculeType()).setDiameter(sigma[i]);
                    sim.potential.setEpsilon(epsilon[i]);
                    sim.potential.setSigma(sigma[i]);
                }
            }
    	    double targetDensity = 0.5;
    	    if(currentSig > 0.5) {
    	        sim.phase.setDensity(targetDensity/Math.pow(currentSig,sim.getSpace().D()));
    	        double size = currentSig*Math.pow(sim.phase.atomCount()/targetDensity,1.0/sim.getSpace().D());
                IVector v = sim.getSpace().makeVector();
                v.E(size);
                sim.phase.getBoundary().setDimensions(v);
                SpaceLattice lattice;
                if (sim.phase.getSpace().D() == 2) {
                    lattice = new LatticeOrthorhombicHexagonal();
                }
                else {
                    lattice = new LatticeCubicFcc();
                }
                Configuration config = new ConfigurationLattice(lattice);
                config.initializeCoordinates(sim.phase);
            }
            if(speciesName.equals("Ideal gas")) sim.getPotentialMaster().setEnabled(sim.potential,false);
            else sim.getPotentialMaster().setEnabled(sim.potential,true);
            simRestart.actionPerformed();
        }
    }
    
    public static class Applet extends javax.swing.JApplet {

	    public void init() {
	        getRootPane().putClientProperty(
	                        "defeatSystemEventQueueCheck", Boolean.TRUE);
		    getContentPane().add(new JouleThomson().panel());
	    }
    }//end of Applet
    
    
    
    public static void main(String[] args) {
        int dim = 3;
        if(args.length > 0) {
            try {
            dim = Integer.parseInt(args[0]);
            } catch(NumberFormatException ex) {
                throw new NumberFormatException(args[0]+" is not a number!");
            }
        }
        
        SimulationGraphic simGraphic = new JouleThomson(Space.getInstance(dim));
		simGraphic.makeAndDisplayFrame("Joule Thomson Experiment");
    }//end of main
        
}