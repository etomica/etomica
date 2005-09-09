package etomica.modules.pistoncylinder;
import java.awt.BorderLayout;
import java.awt.Font;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

import javax.swing.ButtonGroup;
import javax.swing.JComboBox;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JTabbedPane;
import javax.swing.UIManager;

import etomica.action.Action;
import etomica.atom.AtomTypeSphere;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorHistory;
import etomica.data.DataFork;
import etomica.data.DataPump;
import etomica.data.DataSink;
import etomica.data.DataSource;
import etomica.data.DataSourceCountSteps;
import etomica.data.DataSourceScalar;
import etomica.data.DataSourceUniform;
import etomica.data.meter.MeterTemperature;
import etomica.exception.ConfigurationOverlapException;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DeviceBox;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DeviceToggleButton;
import etomica.graphics.DeviceTrioControllerButton;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxesCAE;
import etomica.graphics.DisplayCanvasInterface;
import etomica.graphics.DisplayPhase;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntervalActionAdapter;
import etomica.modifier.Modifier;
import etomica.modifier.ModifierBoolean;
import etomica.modifier.ModifierFunctionWrapper;
import etomica.modifier.ModifierGeneral;
import etomica.potential.P2HardSphere;
import etomica.potential.P2Ideal;
import etomica.potential.P2SquareWell;
import etomica.potential.Potential2HardSphericalWrapper;
import etomica.potential.PotentialGroup;
import etomica.units.Angstrom;
import etomica.units.Bar;
import etomica.units.BaseUnit;
import etomica.units.BaseUnitPseudo3D;
import etomica.units.Dalton;
import etomica.units.Dimension;
import etomica.units.Joule;
import etomica.units.Kelvin;
import etomica.units.Liter;
import etomica.units.Mole;
import etomica.units.Unit;
import etomica.units.UnitRatio;
import etomica.util.Default;


public class PistonCylinderGraphic {
    
    static {
        try {
//            javax.swing.plaf.metal.MetalLookAndFeel.setCurrentTheme(new EtomicaTheme());
//            javax.swing.plaf.metal.MetalLookAndFeel.setCurrentTheme(new etomica.graphics.BlueRoseTheme());
//            javax.swing.UIManager.setLookAndFeel("javax.swing.plaf.metal.MetalLookAndFeel");
//            UIManager.setLookAndFeel("com.sun.java.swing.plaf.motif.MotifLookAndFeel");
            System.out.println(javax.swing.UIManager.getSystemLookAndFeelClassName());
            javax.swing.UIManager.setLookAndFeel(javax.swing.UIManager.getSystemLookAndFeelClassName());
            UIManager.put("TitledBorder.font", new Font("SansSerif", Font.PLAIN, 12));
        } catch(Exception e) {}
    }

    
    public JPanel panel, displayPhasePanel;
    public PistonCylinder pc;
    public Potential2HardSphericalWrapper potentialWrapper;
    public P2HardSphere potentialHS;
    public P2SquareWell potentialSW;
    public P2Ideal potentialIdeal;
    public PotentialGroup potentialGroupHS, potentialGroupSW;
    public DataSourceCountSteps meterCycles;
    public DisplayBox displayCycles; 
    public MeterTemperature thermometer;
    public DisplayPhase displayPhase;
    public DeviceTrioControllerButton controlButtons;
    public ItemListener potentialChooserListener;
    public JComboBox potentialChooser;
    public DeviceSlider scaleSlider, pressureSlider, temperatureSlider;
    public JPanel pressureSliderPanel;
    public MeterPistonDensity densityMeter;
    public DeviceToggleButton fixPistonButton;
    public DisplayPlot plotT, plotD, plotP;
    public Unit tUnit, dUnit, pUnit;
    public final javax.swing.JTabbedPane displayPanel;
    public DeviceBox sigBox, epsBox, lamBox, massBox;
	private DisplayBoxesCAE densityDisplayBox, temperatureDisplayBox, pressureDisplayBox;
    final JRadioButton buttonAdiabatic, buttonIsothermal;
    final JPanel blankPanel = new JPanel();
    public int historyLength;
    public DataSourceWallPressure pressureMeter;
    public int plotUpdateInterval, dataInterval, guiInterval;
    public Unit eUnit;
    public double lambda, epsilon, mass;
    public DeviceSlider doSleepSlider;
    public int repaintSleep = 100;
    public int integratorSleep = 10;
    Default defaults;
    
    public PistonCylinderGraphic() {
        defaults = new Default();
        defaults.blockSize = 100;
        displayPhase = new DisplayPhase(null);
        displayPhase.setColorScheme(new ColorSchemeByType());


        defaults.ignoreOverlap = true;
        defaults.atomSize = 3.0;
        eUnit = new UnitRatio(Joule.UNIT, Mole.UNIT);
        historyLength = 100;
        
        plotUpdateInterval = 10;
        dataInterval = 10;
        guiInterval = 10;

        final int p0 = 500;
        
        lambda = 2.0;
        epsilon = eUnit.toSim(1500.0);
        mass = defaults.atomMass;
        
        //restart action and button
        controlButtons = new DeviceTrioControllerButton();
        
        //adiabatic/isothermal radio button
        ButtonGroup thermalGroup = new ButtonGroup();
        buttonAdiabatic = new JRadioButton("Adiabatic");
        buttonIsothermal = new JRadioButton("Isothermal");
        buttonAdiabatic.setSelected(true);
        thermalGroup.add(buttonAdiabatic);
        thermalGroup.add(buttonIsothermal);
        buttonIsothermal.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent evt) {
                pc.controller.doActionNow( new Action() {
                    public void actionPerformed() {
                        pc.integrator.setIsothermal(buttonIsothermal.isSelected());
                    }
                    public String getLabel() {return "";}
                });
            }
        });
        
        //temperature selector
        temperatureSlider = new DeviceSlider(null);
        temperatureSlider.setShowValues(true);
        temperatureSlider.setEditValues(true);
        temperatureSlider.setMinimum(0);
        temperatureSlider.setMaximum(1000);
        temperatureSlider.getSlider().setMajorTickSpacing(200);
        temperatureSlider.getSlider().setMinorTickSpacing(50);
        temperatureSlider.getSlider().setLabelTable(
        temperatureSlider.getSlider().createStandardLabels(200,100));
        temperatureSlider.setValue(300);

	    //combo box to select potentials
//	    final AtomPairIterator iterator = potentialDisks.iterator();
	    potentialChooser = new javax.swing.JComboBox(new String[] {
	        "Ideal gas", "Repulsion only", "Repulsion and attraction"});
	    potentialChooser.setSelectedIndex(0);

        sigBox = new DeviceBox();
        epsBox = new DeviceBox();
        lamBox = new DeviceBox();
        massBox = new DeviceBox();
        
//        displayPhase.canvas.setDrawBoundary(DisplayCanvasInterface.DRAW_BOUNDARY_NONE);
//        displayPhase.getOriginShift()[0] = thickness;
//        displayPhase.getOriginShift()[1] = -thickness;
	    
	    //slider for scale of display
	    ModifierFunctionWrapper scaleModulator = new ModifierFunctionWrapper(displayPhase, "scale");
	    scaleModulator.setFunction(new etomica.util.Function.Linear(0.01, 0.0));
	    scaleSlider = new DeviceSlider(null, scaleModulator);
	    JPanel scaleSliderPanel = new JPanel(new java.awt.GridLayout(0,1));
	    scaleSliderPanel.add(scaleSlider.graphic());	    
//        scaleSliderPanel.setBorder(new javax.swing.border.TitledBorder("Scale (%)"));
	    scaleSlider.getSlider().addChangeListener(new javax.swing.event.ChangeListener() {
	        public void stateChanged(javax.swing.event.ChangeEvent evt) {
	            if(displayPhase.graphic() != null) displayPhase.graphic().repaint();
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
		
		//add meter and display for current kinetic temperature
		defaults.historyPeriod = 1000;

		thermometer = new MeterTemperature();

        //plot of temperature and density histories
/*		DisplayPlot plotD = new DisplayPlot(this);
		DisplayPlot plotT = new DisplayPlot(this);
		plotD.setDataSources(densityHistory);
        plotD.setUnit(dUnit);
		plotT.setDataSources(temperatureHistory);
		plotT.setUnit(tUnit);
		plotD.setLabel("Density");
		plotT.setLabel("Temperature");
		plotT.getPlot().setYRange(0.0,1500.);
*/		
		//display of averages
/*		DisplayTable table = new DisplayTable(this);
		table.setUpdateInterval(20);
		table.setWhichValues(new MeterAbstract.ValueType[] {
		                MeterAbstract.CURRENT, MeterAbstract.AVERAGE, MeterAbstract.ERROR});
        this.mediator().addMediatorPair(new Mediator.DisplayMeter.NoAction(this.mediator()));
*/				
		
		//pressure device
//        sliderModulator.setFunction(pressureRescale);
//        pressureSlider = new DeviceSelectPressure(controller,integrator);
        pressureSlider = new DeviceSlider(null);
        pressureSlider.setShowValues(true);
        pressureSlider.setEditValues(true);
        pressureSlider.setMinimum(0);
        pressureSlider.setMaximum(1000);
	    pressureSlider.getSlider().setMajorTickSpacing(200);
	    pressureSlider.getSlider().setMinorTickSpacing(50);
	    pressureSlider.getSlider().setLabelTable(
            pressureSlider.getSlider().createStandardLabels(200,100));
	    pressureSlider.setValue(p0);
        
        DeviceSlider repaintSlider = new DeviceSlider(null);
        repaintSlider.setShowValues(true);
        repaintSlider.setEditValues(true);
        repaintSlider.setMinimum(0);
        repaintSlider.setMaximum(1000);
        repaintSlider.getSlider().setMajorTickSpacing(200);
        repaintSlider.getSlider().setLabelTable(
                repaintSlider.getSlider().createStandardLabels(200,100));
        repaintSlider.setValue(repaintSleep);
        repaintSlider.setModifier(new Modifier() {
            public String getLabel() {return "";}
            public Dimension getDimension() {return Dimension.UNDEFINED;}
            public void setValue(double v) {
                repaintSleep = (int)v;
            }
            public double getValue() {return (double)repaintSleep;}
        });
        
        doSleepSlider = new DeviceSlider(null);
        doSleepSlider.setShowValues(true);
        doSleepSlider.setEditValues(true);
        doSleepSlider.setMinimum(0);
        doSleepSlider.setMaximum(100);
        doSleepSlider.getSlider().setMajorTickSpacing(20);
        doSleepSlider.getSlider().setLabelTable(
                doSleepSlider.getSlider().createStandardLabels(20,10));
        doSleepSlider.setValue(integratorSleep);

        //set-pressure history
//        etomica.MeterScalar pressureSetting = new MeterDatumSourceWrapper(pressureSlider.getModulator());
//        pressureSetting.setHistorying(true);
//        History pHistory = pressureSetting.getHistory();
//        pHistory.setLabel("Set Pressure");
/*        DisplayPlot plotP = new DisplayPlot(this);
        plotP.setDataSources(pHistory);
        plotP.setUnit(pUnit);
        plotP.setLabel("Pressure"); 
        plotP.getPlot().setYRange(0.0, 1500.);
*/        
        //measured pressure on piston
        //wrap it in a MeterDatumSource wrapper because we want to average pressure
        //over a longer interval than is used by other meters.  By wrapping it
        //we can still have history synchronized with others
//        Atom piston = ((SpeciesAgent)pistonCylinder.getAgent(phase)).node.firstLeafAtom();
//        piston.coord.setMass(pistonMass);
//        etomica.MeterScalar pressureMeter = ((AtomType.Wall)piston.type).new MeterPressure(this);
//        pressureMeter.setUpdateInterval(10);
//        pressureMeter.setFunction(pressureScale);
//        etomica.MeterDatumSourceWrapper pressureMeterWrapper = new MeterDatumSourceWrapper(pressureMeter);
////        pressureMeterWrapper.setFunction(pressureRescale);
//        pressureMeterWrapper.setWhichValue(MeterAbstract.MOST_RECENT);
//        pressureMeterWrapper.setHistorying(true);
//        History pMeterHistory = pressureMeterWrapper.getHistory();
//        pMeterHistory.setLabel("Pressure ("+pUnit.symbol()+")");
        
//        DisplayPlot plot = new DisplayPlot(this);
//        plot.setDataSources(new DataSource[] {
//                densityHistory, temperatureHistory, pMeterHistory, pHistory});
//        plot.setLabel("History");
//        plot.getPlot().setYLabel("");
//        plot.getPlot().setYRange(0.0, 1500.);
//        plot.setYUnit(new Unit[] {dadUnit, tUnit, pUnit, pUnit});
        
        fixPistonButton = new DeviceToggleButton(null);

        plotD = new DisplayPlot();
        plotT = new DisplayPlot();
        plotP = new DisplayPlot();


        //************* Lay out components ****************//
        
        panel = new JPanel();
        panel.setLayout(new java.awt.BorderLayout());      

        //tabbed pane for the big displays
    	displayPanel = new javax.swing.JTabbedPane();
    	displayPhasePanel = new javax.swing.JPanel(new java.awt.BorderLayout());
    	displayPhasePanel.add(scaleSliderPanel,java.awt.BorderLayout.EAST);
        
        JPanel plotPanel = new JPanel(new java.awt.GridLayout(0,1));
        plotPanel.add(plotD.graphic());
        plotPanel.add(plotT.graphic());
        plotPanel.add(plotP.graphic());
        displayPanel.add("Plots",new javax.swing.JScrollPane(plotPanel));

        JPanel startPanel = (JPanel)controlButtons.graphic();
        java.awt.GridBagConstraints gbc0 = new java.awt.GridBagConstraints();
        startPanel.setBorder(new javax.swing.border.TitledBorder("Control"));
        gbc0.gridx = 0; gbc0.gridy = 0;
        gbc0.gridx = 0; gbc0.gridy = 2; gbc0.gridwidth = 2;
        startPanel.add(fixPistonButton.graphic(null), gbc0);
        startPanel.setLayout(new GridLayout(2,2));

        //panel for the temperature control/display
        JPanel temperaturePanel = new JPanel(new java.awt.GridBagLayout());
        temperaturePanel.setBorder(new javax.swing.border.TitledBorder("Set Temperature (K)"));
        java.awt.GridBagConstraints gbc1 = new java.awt.GridBagConstraints();
        gbc1.gridx = 0;  gbc1.gridy = 1;
        gbc1.gridwidth = 1;
        temperaturePanel.add(buttonAdiabatic,gbc1);
        gbc1.gridx = 1;  gbc1.gridy = 1;
        gbc1.gridwidth = 1;
        temperaturePanel.add(buttonIsothermal,gbc1);
        gbc1.gridx = 0;  gbc1.gridy = 2;
        gbc1.gridwidth = 2;
        temperaturePanel.add(temperatureSlider.graphic(),gbc1);
        
        //panel for pressure slider
        pressureSliderPanel = new JPanel(new java.awt.GridLayout(0,1));
        pressureSlider.setShowBorder(false);
        pressureSliderPanel.add(pressureSlider.graphic());
        
        //panel for all the controls
        JPanel dimensionPanel = new JPanel(new GridLayout(1,0));
        ButtonGroup dimensionGroup = new ButtonGroup();
        final JRadioButton button2D = new JRadioButton("2D");
        JRadioButton button3D = new JRadioButton("3D");
        button2D.setSelected(true);
        dimensionGroup.add(button2D);
        dimensionGroup.add(button3D);
        dimensionPanel.add(button2D);
        dimensionPanel.add(button3D);
        button2D.addItemListener(new ItemListener() {
           public void itemStateChanged(ItemEvent evt) {
               if(button2D.isSelected()) {
                   setSimulation(new PistonCylinder(2, defaults));
               } else {
                   setSimulation(new PistonCylinder(3, defaults));
               }
           }
        });

        JPanel controlPanel = new JPanel(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gbc2 = new java.awt.GridBagConstraints();
        gbc2.gridx = 0;
        gbc2.gridy = java.awt.GridBagConstraints.RELATIVE;
        controlPanel.add(dimensionPanel);
        controlPanel.add(startPanel,gbc2);
        java.awt.GridBagConstraints gbc3 = new java.awt.GridBagConstraints();
        gbc3.gridx = java.awt.GridBagConstraints.RELATIVE;
        gbc3.gridy = 0;
        displayCycles = new DisplayBox();
        controlPanel.add(displayCycles.graphic(),gbc2);

        JPanel repaintSliderPanel = new JPanel(new java.awt.GridLayout(0,1));
        repaintSlider.setShowBorder(false);
        repaintSliderPanel.setBorder(new javax.swing.border.TitledBorder("Repaint delay"));
        repaintSliderPanel.add(repaintSlider.graphic());

        JPanel doSleepSliderPanel = new JPanel(new java.awt.GridLayout(0,1));
        repaintSlider.setShowBorder(false);
        doSleepSliderPanel.setBorder(new javax.swing.border.TitledBorder("Integrator step delay"));
        doSleepSliderPanel.add(doSleepSlider.graphic());
        
        JPanel statePanel = new JPanel(new GridBagLayout());
        statePanel.add(temperaturePanel, gbc2);
        statePanel.add(pressureSliderPanel, gbc2);

        JPanel potentialPanel = new JPanel(new GridBagLayout());
        potentialPanel.add(potentialChooser,gbc2);
	    potentialPanel.setBorder(new javax.swing.border.TitledBorder("Potential selection"));
	    JPanel parameterPanel = new JPanel(new GridLayout(0,1));
        parameterPanel.add(sigBox.graphic());
        parameterPanel.add(epsBox.graphic());
        parameterPanel.add(lamBox.graphic());
        parameterPanel.add(massBox.graphic());
        potentialPanel.add(parameterPanel,gbc2);
        
        JPanel guiPanel = new JPanel(new GridBagLayout());
        guiPanel.add(repaintSliderPanel, gbc2);
        guiPanel.add(doSleepSliderPanel, gbc2);
        
        JTabbedPane setupPanel = new JTabbedPane();
        setupPanel.add(statePanel, "State");
        setupPanel.add(potentialPanel, "Potential");
        setupPanel.add(guiPanel, "Controls");
        
        //panel for the density, temperature, and pressure displays
        JPanel dataPanel = new JPanel(new GridBagLayout());
        
        densityDisplayBox = new DisplayBoxesCAE();
        densityDisplayBox.setLabel("Density (mol/L)");
        densityDisplayBox.setLabelType(DisplayBox.BORDER);
        dataPanel.add(densityDisplayBox.graphic(),gbc2);
        
        temperatureDisplayBox = new DisplayBoxesCAE();
        temperatureDisplayBox.setLabel("Temperature (K)");
        temperatureDisplayBox.setLabelType(DisplayBox.BORDER);
        dataPanel.add(temperatureDisplayBox.graphic(),gbc2);
        
        pressureDisplayBox = new DisplayBoxesCAE();
        pressureDisplayBox.setLabel("Pressure (bar)");
        pressureDisplayBox.setLabelType(DisplayBox.BORDER);
        dataPanel.add(pressureDisplayBox.graphic(),gbc2);
        
        JPanel leftPanel = new JPanel(new GridBagLayout());
        
        leftPanel.add(controlPanel, gbc2);
        leftPanel.add(setupPanel, gbc2);
        leftPanel.add(dataPanel, gbc2);
        
        panel.add(leftPanel, BorderLayout.WEST);
        panel.add(displayPanel, BorderLayout.EAST);

        Thread repainter = new Thread() {
            public void run() {
                while (true) {
                    panel.repaint();
                    try{Thread.sleep(repaintSleep);}
                    catch(InterruptedException e){}
                }
            }
        };
        repainter.start();
        
        setSimulation(new PistonCylinder(2, defaults));
    }
    
    public void setSimulation(PistonCylinder sim) {
        boolean pistonHeld = true;
        if (pc != null) {
            pistonHeld = pc.pistonPotential.isStationary();
            lambda = potentialSW.getLambda();
            epsilon = potentialSW.getEpsilon();
            mass = pc.species.getMass();
            pc.getController().halt();
        }
        pc = sim;
        pc.species.setMass(mass);
        int D = pc.space.D();


        tUnit = Kelvin.UNIT;

        if (pc.space.D() == 2) {
            dUnit = new UnitRatio(Mole.UNIT, 
                    new BaseUnitPseudo3D.Volume(Liter.UNIT));
            pUnit = new BaseUnitPseudo3D.Pressure(Bar.UNIT);
        }
        else {
            dUnit = new UnitRatio(Mole.UNIT, Liter.UNIT);
            pUnit = Bar.UNIT;
        }

        // set up GUI
        BaseUnit.Length.Sim.TO_PIXELS = 600/pc.phase.boundary().dimensions().x(1);
        pc.ai.setDoSleep(true);
        pc.integrator.setTimeStep(0.5);
        pc.integrator.removeAllListeners();

        pc.wallPotential.setLongWall(0,true,true);  // left wall
        pc.wallPotential.setLongWall(0,false,true); // right wall
        // skip top wall
        pc.wallPotential.setLongWall(1,false,false);// bottom wall
        pc.wallPotential.setPhase(pc.phase);  // so it has a boundary
        
        if (displayPhase.graphic() != null) {
            displayPhasePanel.remove(displayPhase.graphic());
            displayPanel.remove(displayPhasePanel);
            displayPanel.remove(blankPanel);
        }
        if (D == 2) {
            pc.ai.setSleepPeriod(integratorSleep);
            displayPanel.insertTab(displayPhase.getLabel(), null, displayPhasePanel, "", 0);
//            displayPanel.add(displayPhase.getLabel(), displayPhasePanel);
            displayPhase.setPhase(pc.phase);
            displayPhase.setAlign(1,DisplayPhase.BOTTOM);
            displayPhase.canvas.setDrawBoundary(DisplayCanvasInterface.DRAW_BOUNDARY_NONE);
            displayPhase.getDrawables().clear();
            displayPhase.addDrawable(pc.pistonPotential);
            displayPhase.addDrawable(pc.wallPotential);
            displayPhasePanel.add(displayPhase.graphic(),java.awt.BorderLayout.WEST);
//            pc.integrator.addIntervalListener(displayPhase);
        } else {
            pc.ai.setSleepPeriod(integratorSleep);
            displayPanel.add("Run Faster", blankPanel);
        }
        scaleSlider.setController(pc.controller);

        //  control panel
        ModifierBoolean fixPistonModulator = new ModifierBoolean() {
            public void setBoolean(boolean b) {
                pc.pistonPotential.setStationary(b);
            }
            public boolean getBoolean() {
                return pc.pistonPotential.isStationary();
            }
        };
        fixPistonButton.setController(pc.controller);
        fixPistonButton.setModifier(fixPistonModulator, "Release piston", "Hold piston");
        fixPistonButton.setPostAction(new ActionPistonUpdate(pc.integrator));
        fixPistonButton.setState(pistonHeld);

        meterCycles = new DataSourceCountSteps();
        pc.integrator.addListener(meterCycles);
        displayCycles.setPrecision(6);
        DataPump pump= new DataPump(meterCycles,displayCycles);
        new IntervalActionAdapter(pump,pc.integrator);
        displayCycles.setLabel("Integrator steps");
        controlButtons.setSimulation(pc);
        controlButtons.getRestartButton().setPreAction(new Action() {
            public String getLabel() {
                return "Piston Reset";
            }
            public void actionPerformed() {
                pc.pistonPotential.setWallPosition(0.0);
            }
        });

        
        //  state panel
        pc.integrator.setIsothermal(buttonIsothermal.isSelected());
        pc.integrator.setTemperature(tUnit.toSim(temperatureSlider.getValue()));
        temperatureSlider.setUnit(tUnit);
        temperatureSlider.setModifier(new ModifierGeneral(pc.integrator,"temperature"));
        temperatureSlider.setController(pc.getController());

        potentialSW = new P2SquareWell(pc.space,defaults.atomSize,lambda,epsilon,defaults.ignoreOverlap);
        potentialHS = new P2HardSphere(pc.space,defaults.atomSize,defaults.ignoreOverlap);
        potentialIdeal = new P2Ideal(pc.space);
        
        if(potentialChooserListener != null) potentialChooser.removeItemListener(potentialChooserListener);
        
        potentialChooserListener = new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent evt) {
                if(evt.getStateChange() == java.awt.event.ItemEvent.DESELECTED) return;
                setPotential((String)evt.getItem());
            }
        };
        potentialChooser.addItemListener(potentialChooserListener);
        setPotential((String)potentialChooser.getSelectedItem());

        ModifierAtomDiameter sigModifier = new ModifierAtomDiameter();
        ModifierGeneral epsModifier = new ModifierGeneral(potentialSW, "epsilon");
        ModifierGeneral lamModifier = new ModifierGeneral(potentialSW, "lambda");
        ModifierGeneral massModifier = new ModifierGeneral(pc.species,"mass");
        sigBox.setModifier(sigModifier);
        sigBox.setLabel("Core Diameter ("+Angstrom.UNIT.symbol()+")");
        epsBox.setUnit(eUnit);
        epsBox.setModifier(epsModifier);
        lamBox.setModifier(lamModifier);
        massBox.setModifier(massModifier);
        massBox.setUnit(Dalton.UNIT);
        sigBox.setController(pc.getController());
        epsBox.setController(pc.getController());
        lamBox.setController(pc.getController());
        massBox.setController(pc.getController());
        
        pressureSlider.setUnit(pUnit);
        pressureSliderPanel.setBorder(new javax.swing.border.TitledBorder("Set Pressure ("+pUnit.toString()+")"));
        Dimension pDim = (D==2) ? Dimension.PRESSURE2D : Dimension.PRESSURE;
        pc.pistonPotential.setPressure(pUnit.toSim(pressureSlider.getValue()));
        pressureSlider.setModifier(new ModifierPistonPressure(pc.pistonPotential,pDim));
        pressureSlider.setPostAction(new ActionPistonUpdate(pc.integrator));
        pressureSlider.setController(pc.getController());

        doSleepSlider.setModifier(new Modifier() {
            public String getLabel() {return "";}
            public Dimension getDimension() {return Dimension.UNDEFINED;}
            public void setValue(double v) {
                integratorSleep = (int)v;
                if (integratorSleep > 0) {
                    pc.ai.setDoSleep(true);
                    pc.ai.setSleepPeriod(integratorSleep);
                }
                else {
                    pc.ai.setDoSleep(false);
                }
            }
            public double getValue() {return integratorSleep;}
        });
        
        //  data panel
        // the data channel sending the DataTable should be cut off (and hopefully garbage collected)
        plotD.getDataTable().reset();
        plotT.getDataTable().reset();
        plotP.getDataTable().reset();
        
        thermometer.setPhase(pc.phase);
        AccumulatorHistory temperatureHistory = new AccumulatorHistory();
        temperatureHistory.setHistoryLength(historyLength);
        AccumulatorAverage temperatureAvg = new AccumulatorAverage(sim);
        pump = new DataPump(thermometer,new DataFork(new DataSink[]{temperatureHistory,temperatureAvg}));
        IntervalActionAdapter adapter = new IntervalActionAdapter(pump,pc.integrator);
        adapter.setActionInterval(dataInterval);
        temperatureHistory.addDataSink(plotT.getDataTable());
        temperatureDisplayBox.setAccumulator(temperatureAvg);
        temperatureDisplayBox.setUnit(tUnit);
        
        DataSource targetTemperatureDataSource = new DataSourceScalar("Temperature",Dimension.TEMPERATURE) {
            public double getDataAsScalar() {
                return temperatureSlider.getModifier().getValue();
            }
        };
        AccumulatorHistory targetTemperatureHistory = new AccumulatorHistory();
        targetTemperatureHistory.setHistoryLength(historyLength);
        DataPump targetTemperatureDataPump = new DataPump(targetTemperatureDataSource, targetTemperatureHistory);
        adapter = new IntervalActionAdapter(targetTemperatureDataPump,pc.integrator);
        adapter.setActionInterval(dataInterval);
        targetTemperatureHistory.addDataSink(plotT.getDataTable());

        pressureMeter = new DataSourceWallPressure(pc.space,pc.pistonPotential,pc.integrator);
        pressureMeter.setPhase(pc.phase);
        AccumulatorHistory pressureHistory = new AccumulatorHistory();
        pressureHistory.setHistoryLength(historyLength);
        AccumulatorAverage pressureAvg = new AccumulatorAverage(sim);
        pump = new DataPump(pressureMeter, new DataFork(new DataSink[]{pressureHistory,pressureAvg}));
        adapter = new IntervalActionAdapter(pump,pc.integrator);
        adapter.setActionInterval(dataInterval);
        pressureHistory.addDataSink(plotP.getDataTable());
        pressureDisplayBox.setAccumulator(pressureAvg);
        pressureDisplayBox.setUnit(pUnit);

        DataSource targetPressureDataSource = new DataSourceScalar("Pressure",Dimension.PRESSURE) {
            public double getDataAsScalar() {
                return pUnit.toSim(pressureSlider.getValue());
            }
        };
        AccumulatorHistory targetPressureHistory = new AccumulatorHistory();
        targetPressureHistory.setHistoryLength(historyLength);
        pump = new DataPump(targetPressureDataSource, targetPressureHistory);
        adapter = new IntervalActionAdapter(pump,pc.integrator);
        adapter.setActionInterval(dataInterval);
        targetPressureHistory.addDataSink(plotP.getDataTable());

        densityMeter = new MeterPistonDensity(pc.pistonPotential,1,defaults.atomSize);
        densityMeter.setPhase(pc.phase);
        AccumulatorHistory densityHistory = new AccumulatorHistory();
        densityHistory.setHistoryLength(historyLength);
        AccumulatorAverage densityAvg = new AccumulatorAverage(sim);
        pump = new DataPump(densityMeter,new DataFork(new DataSink[]{densityAvg, densityHistory}));
        adapter = new IntervalActionAdapter(pump,pc.integrator);
        adapter.setActionInterval(dataInterval);
        densityHistory.addDataSink(plotD.getDataTable());
        densityDisplayBox.setAccumulator(densityAvg);
        densityDisplayBox.setUnit(dUnit);
        
        plotD.setUnit(dUnit);
        plotT.setUnit(tUnit);
        plotP.setUnit(pUnit);

        plotP.getPlot().setTitle("Pressure ("+pUnit.symbol()+")");
        plotT.getPlot().setTitle("Temperature ("+tUnit.symbol()+")");
        plotD.getPlot().setTitle("Density ("+dUnit.symbol()+")");
        
        plotT.getPlot().setXRange(0, historyLength);
        plotP.getPlot().setXRange(0, historyLength);
        plotD.getPlot().setXRange(0, historyLength);
        plotT.getPlot().setYRange(0, 1000.);
        plotP.getPlot().setYRange(0, 1000.);
        DataSourceUniform xSource = new DataSourceUniform("History",Dimension.QUANTITY,historyLength, 1, historyLength);
        plotT.setXSource(xSource);
        plotP.setXSource(xSource);
        plotD.setXSource(xSource);

        java.awt.Dimension d = plotT.getPlot().getPreferredSize();
        d.height = 230;
        plotT.getPlot().setSize(d);
        plotP.getPlot().setSize(d);
        plotD.getPlot().setSize(d);
    }
    
    public void setPotential(String potentialDesc) {
        final boolean HS = potentialDesc.equals("Repulsion only"); 
        final boolean SW = potentialDesc.equals("Repulsion and attraction"); 
        pc.controller.doActionNow( new Action() {
            public void actionPerformed() {
                if (HS) {
                    pc.potentialWrapper.setPotential(potentialHS);
                }
                else if (SW) {
                    pc.potentialWrapper.setPotential(potentialSW);
                }
                else {
                    pc.potentialWrapper.setPotential(potentialIdeal);
                }
                try {
                    pc.integrator.reset();
                } catch(ConfigurationOverlapException e) {}
            }
            public String getLabel() {return "";}
        });
    }
    
    private class ModifierAtomDiameter implements Modifier {

        public void setValue(double d) {
            defaults.atomSize = d;
            //assume one type of atom
            ((AtomTypeSphere)pc.phase.firstAtom().type).setDiameter(d);
            PistonCylinderGraphic.this.densityMeter.setAtomDiameter(d);
            PistonCylinderGraphic.this.potentialHS.setCollisionDiameter(d);
            PistonCylinderGraphic.this.potentialSW.setCoreDiameter(d);
            pc.pistonPotential.setCollisionRadius(0.5*d);
            pc.wallPotential.setCollisionRadius(0.5*d);
            displayPhase.repaint();
        }

        public double getValue() {
            return defaults.atomSize;
        }

        public Dimension getDimension() {
            return Dimension.LENGTH;
        }
        
        public String getLabel() {
            return "Atom Diameter";
        }
        
        public String toString() {
            return getLabel();
        }
    }
    
    public static void main(String[] args) {
        PistonCylinderGraphic sim = new PistonCylinderGraphic();
		SimulationGraphic.makeAndDisplayFrame(sim.panel);
//		sim.phase.reset();
 //       sim.controller1.start();
    }//end of main
    
    public static class Applet extends javax.swing.JApplet {

	    public void init() {
		    getContentPane().add(new PistonCylinderGraphic().panel);
	    }
    }//end of Applet
}//end of PistonCylinderGraphic class


