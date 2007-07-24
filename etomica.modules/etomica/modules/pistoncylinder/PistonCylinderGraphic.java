package etomica.modules.pistoncylinder;
import java.awt.BorderLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.event.ItemListener;
import java.util.ArrayList;

import javax.swing.JComboBox;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.border.TitledBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import etomica.action.Action;
import etomica.action.IntegratorReset;
import etomica.action.SimulationRestart;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.AtomTypeSphere;
import etomica.chem.elements.ElementSimple;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.AccumulatorHistory;
import etomica.data.DataFork;
import etomica.data.DataPump;
import etomica.data.DataSink;
import etomica.data.DataSource;
import etomica.data.DataSourceCountTime;
import etomica.data.DataSourceScalar;
import etomica.data.DataTag;
import etomica.data.meter.MeterDensity;
import etomica.data.meter.MeterTemperature;
import etomica.exception.ConfigurationOverlapException;
import etomica.graphics.ActionConfigWindow;
import etomica.graphics.ActionVelocityWindow;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DeviceBox;
import etomica.graphics.DeviceButton;
import etomica.graphics.DeviceNSelector;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DeviceThermoSlider;
import etomica.graphics.DeviceToggleButton;
import etomica.graphics.DeviceTrioControllerButton;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayCanvasInterface;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.integrator.IIntegrator;
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
import etomica.units.CompoundUnit;
import etomica.units.Dalton;
import etomica.units.Dimension;
import etomica.units.Joule;
import etomica.units.Kelvin;
import etomica.units.Length;
import etomica.units.Liter;
import etomica.units.Meter;
import etomica.units.Mole;
import etomica.units.Pixel;
import etomica.units.Prefix;
import etomica.units.PrefixedUnit;
import etomica.units.Pressure;
import etomica.units.Temperature;
import etomica.units.Time;
import etomica.units.Undefined;
import etomica.units.Unit;
import etomica.units.UnitRatio;
import etomica.units.systems.MKS;


public class PistonCylinderGraphic extends SimulationPanel {
    
	private static final String APP_NAME = "Piston Cylinder";
    public JPanel displayBoxPanel;
    public PistonCylinder pc;
    public Potential2HardSphericalWrapper potentialWrapper;
    public P2HardSphere potentialHS;
    public P2SquareWell potentialSW;
    public P2Ideal potentialIdeal;
    public PotentialGroup potentialGroupHS, potentialGroupSW;
    public DataSourceCountTime meterCycles;
    public DisplayTextBox displayCycles;
    public MeterTemperature thermometer;
    public DisplayBox displayBox;
    public DeviceTrioControllerButton controlButtons;
    public DeviceButton configButton, velocityButton;
    public ItemListener potentialChooserListener;
    public JComboBox potentialChooser;
    DeviceThermoSlider tempSlider;
    public DeviceSlider scaleSlider, pressureSlider;
    public DeviceNSelector nSlider;
    public JPanel pressureSliderPanel;
    public MeterDensity densityMeter;
    public DeviceToggleButton fixPistonButton;
    public DisplayPlot plotT, plotD, plotP;
    public DisplayPlot plotRDF;
    public Unit tUnit, dUnit, pUnit;
    public DeviceBox sigBox, epsBox, lamBox, massBox;
	private DisplayTextBoxesCAE densityDisplayTextBox, temperatureDisplayTextBox, pressureDisplayTextBox;
    public JPanel blankPanel = new JPanel();
    public JScrollPane plotsPane;
    public int historyLength;
    public DataSourceWallPressure pressureMeter;
    public int dataInterval;
    public Unit eUnit;
    public double lambda, epsilon, mass, sigma;
    public DeviceSlider doSleepSlider, integratorTimeStepSlider;
    public int repaintSleep = 100;
    public int integratorSleep = 10;
    public int integratorSleep3D = 0;
    private boolean initialized;
    private boolean doConfigButton = false;
    private boolean doRDF = false;
    private boolean showTimeControls = false;
    private JPanel guiPanel;

    public PistonCylinderGraphic() {
    	super(APP_NAME);
    }

    /**
     * Enable/disable button to show coordinates.  Must be called before init.
     */
    public void setDoConfigButton(boolean newDoConfigButton) {
        if (initialized) {
            throw new RuntimeException("Already initialized");
        }
        doConfigButton = newDoConfigButton;
    }
    
    /**
     * Enable/disable RDF tab.  Must be called before init.
     */
    public void setDoRDF(boolean newDoRDF) {
        if (initialized) {
            throw new RuntimeException("Already initialized");
        }
        doRDF = newDoRDF;
    }
    
    public void init() {
        initialized = true;

        pc = new PistonCylinder(2);

        displayBox = new DisplayBox(null);
        displayBox.setColorScheme(new ColorSchemeByType());

        eUnit = new UnitRatio(Joule.UNIT, Mole.UNIT);
        historyLength = 100;
        
        dataInterval = 10;

        final int p0 = 500;
        
        lambda = 2.0;
        epsilon = eUnit.toSim(1500.0);
        mass = 40;
        sigma = 4.0;

        GridBagConstraints vertGBC = SimulationPanel.getVertGBC();
        GridBagConstraints horizGBC = SimulationPanel.getHorizGBC();

        // Dimension Selection
        // Currently commented out so that user can only run a 2D simulation
/*
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
                   setSimulation(new PistonCylinder(2));
            } else {
                 setSimulation(new PistonCylinder(3));
             }
         }
      });
*/

        // Control buttons
        controlButtons = new DeviceTrioControllerButton();
        fixPistonButton = new DeviceToggleButton(null);

        JPanel startPanel = (JPanel)controlButtons.graphic();
        GridBagConstraints gbc0 = new GridBagConstraints();
        startPanel.setBorder(new TitledBorder("Control"));
        gbc0.gridx = 0; gbc0.gridy = 0;
        gbc0.gridx = 0; gbc0.gridy = 2; gbc0.gridwidth = 2;
        startPanel.add(fixPistonButton.graphic(null), gbc0);
        startPanel.setLayout(new GridLayout(2,2));

        // Simulation Time
        displayCycles = new DisplayTextBox();

        // Add dimension buttons panel, controls panel,
        // simulation timebox and config buttons, if requested,
        // onto a panel
        JPanel controlsPanel = new JPanel(new GridBagLayout());
//      controlsPanel.add(dimensionPanel);
        controlsPanel.add(startPanel,vertGBC);

        // Show buttons (added to panel 
        if (doConfigButton) {
            JPanel configPanel = new JPanel(new GridBagLayout());
            controlsPanel.add(configPanel,vertGBC);

            configButton = new DeviceButton(null);
            configButton.setLabel("Show Config");
            configPanel.add(configButton.graphic(),horizGBC);

            velocityButton = new DeviceButton(null);
            velocityButton.setLabel("Show Velocities");
            configPanel.add(velocityButton.graphic(),horizGBC);
        }

        //
        // State tabbed pane page
        //

        tempSlider = new DeviceThermoSlider(pc.controller);
        tempSlider.setShowValues(true);
        tempSlider.setEditValues(true);
        tempSlider.setMinimum(0);
        tempSlider.setMaximum(1000);
        tempSlider.setAdiabatic();
        tempSlider.setSliderMajorValues(4);
        tempSlider.setTemperature(300);
        tempSlider.setController(pc.getController());
        tempSlider.setIntegrator(pc.integrator);

		//pressure device
        pressureSlider = new DeviceSlider(null);
        pressureSlider.setShowValues(true);
        pressureSlider.setEditValues(true);
        pressureSlider.setMinimum(0);
        pressureSlider.setMaximum(1000);
        pressureSlider.setNMajor(4);
	    pressureSlider.setValue(p0);

        // panel for pressure control / display
        pressureSliderPanel = new JPanel(new java.awt.GridLayout(0,1));
        pressureSlider.setShowBorder(false);
        pressureSliderPanel.add(pressureSlider.graphic());

        JPanel nSliderPanel = null;


        nSlider = new DeviceNSelector();
        nSlider.setLabel("Number of atoms");
        nSlider.setShowBorder(true);
        nSlider.setEditValues(true);
        nSlider.setShowValues(true);
        // add a listener to adjust the thermostat interval for different
        // system sizes (since we're using ANDERSEN_SINGLE).  Smaller systems 
        // don't need as much thermostating.
        nSlider.getSlider().addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent evt) {
                int n = (int)nSlider.getValue();
                if(n == 0) {
                    pc.integrator.setThermostatInterval(200);
                }
                else {
                	pc.integrator.setThermostatInterval(200/n);
                }
            }
        });

        nSliderPanel = new JPanel(new GridLayout(0,1));
        nSliderPanel.setBorder(new TitledBorder("Number of Molecules"));
        nSlider.setShowBorder(false);
        nSlider.setNMajor(4);
        nSliderPanel.add(nSlider.graphic());


        // Add all state page sub panels onto a single panel
        JPanel statePanel = new JPanel(new GridBagLayout());
        GridBagConstraints gbc2 = new GridBagConstraints();
        gbc2.gridx = 0;  gbc2.gridy = 0;
        statePanel.add(tempSlider.graphic(), gbc2);
        gbc2.gridx = 0;  gbc2.gridy = 1;
        statePanel.add(pressureSliderPanel, gbc2);

        gbc2.gridx = 0;  gbc2.gridy = 2;
        statePanel.add(nSliderPanel, gbc2);


        //
        // Potential tabbed pane page
        //
        
	    //combo box to select potentials
	    potentialChooser = new javax.swing.JComboBox(new String[] {
	        "Ideal gas", "Repulsion only", "Repulsion and attraction"});
	    potentialChooser.setSelectedIndex(0);

        sigBox = new DeviceBox();
        epsBox = new DeviceBox();
        lamBox = new DeviceBox();
        massBox = new DeviceBox();

        JPanel potentialPanel = new JPanel(new GridBagLayout());
        potentialPanel.add(potentialChooser,vertGBC);
	    potentialPanel.setBorder(new TitledBorder("Potential selection"));
	    JPanel parameterPanel = new JPanel(new GridLayout(0,1));
        parameterPanel.add(sigBox.graphic());
        parameterPanel.add(epsBox.graphic());
        parameterPanel.add(lamBox.graphic());
        parameterPanel.add(massBox.graphic());
        potentialPanel.add(parameterPanel,vertGBC);
        
        //
        // Controls tabbed pane page
        //

        if(showTimeControls == true) {
	        // Repaint delay slider
	        DeviceSlider repaintSlider = new DeviceSlider(null);
	        //XXX ugh, see bug 49
	        repaintSlider.setUnit(Time.SIM_UNIT);
	        repaintSlider.setShowValues(false);
	        repaintSlider.setEditValues(true);
	        repaintSlider.setMinimum(0);
	        repaintSlider.setMaximum(1000);
	        repaintSlider.setNMajor(4);
	        repaintSlider.setValue(repaintSleep);
	        repaintSlider.setModifier(new Modifier() {
	            public String getLabel() {return "";}
	            public Dimension getDimension() {return Undefined.DIMENSION;}
	            public void setValue(double v) {
	                repaintSleep = (int)v;
	            }
	            public double getValue() {return repaintSleep;}
	        });
	        JPanel repaintSliderPanel = new JPanel(new GridLayout(0,1));
	        repaintSlider.setShowBorder(false);
	        repaintSliderPanel.setBorder(new TitledBorder("Repaint delay (ms)"));
	        repaintSliderPanel.add(repaintSlider.graphic());
	        
	        // Integrator step delay slider
	        doSleepSlider = new DeviceSlider(null);
	        //XXX ugh, see bug 49
	        doSleepSlider.setUnit(Time.SIM_UNIT);
	        doSleepSlider.setShowValues(false);
	        doSleepSlider.setEditValues(true);
	        doSleepSlider.setMinimum(0);
	        doSleepSlider.setMaximum(100);
	        doSleepSlider.setNMajor(5);
	        JPanel doSleepSliderPanel = new JPanel(new GridLayout(0,1));
	        doSleepSlider.setShowBorder(false);
	        doSleepSliderPanel.setBorder(new TitledBorder("Integrator step delay (ms)"));
	        doSleepSliderPanel.add(doSleepSlider.graphic());
	        
	        // Integrator time step slider
	        integratorTimeStepSlider = new DeviceSlider(null);
	        integratorTimeStepSlider.setShowValues(false);
	        integratorTimeStepSlider.setEditValues(true);
	        integratorTimeStepSlider.setPrecision(2);
	        integratorTimeStepSlider.setMinimum(0.0);
	        integratorTimeStepSlider.setMaximum(5);
	        integratorTimeStepSlider.setNMajor(5);
	        integratorTimeStepSlider.setValue(0.5);
	        JPanel integratorTimeStepSliderPanel = new JPanel(new GridLayout(0,1));
	        repaintSlider.setShowBorder(false);
	        integratorTimeStepSliderPanel.setBorder(new TitledBorder("Integrator time step (ps)"));
	        integratorTimeStepSliderPanel.add(integratorTimeStepSlider.graphic());
	
	        guiPanel = new JPanel(new GridBagLayout());
	        guiPanel.add(repaintSliderPanel, vertGBC);
	        guiPanel.add(doSleepSliderPanel, vertGBC);
	        guiPanel.add(integratorTimeStepSliderPanel, vertGBC);
        }

        //
        // Tabbed pane for state, potential, controls pages
        //
        JTabbedPane setupPanel = new JTabbedPane();
        setupPanel.add(statePanel, "State");
        setupPanel.add(potentialPanel, "Potential");
        if(showTimeControls == true) {
            setupPanel.add(guiPanel, "Controls");
        }

        // Density value display
        densityDisplayTextBox = new DisplayTextBoxesCAE();
        densityDisplayTextBox.setLabelType(DisplayTextBox.LabelType.BORDER);

        // Temperature value display
        temperatureDisplayTextBox = new DisplayTextBoxesCAE();
        temperatureDisplayTextBox.setLabel("Temperature (K)");
        temperatureDisplayTextBox.setLabelType(DisplayTextBox.LabelType.BORDER);

        // Pressure value display
        pressureDisplayTextBox = new DisplayTextBoxesCAE();
        pressureDisplayTextBox.setLabelType(DisplayTextBox.LabelType.BORDER);

        //panel for the density, temperature, and pressure displays
        JPanel dataPanel = new JPanel(new GridBagLayout());
        dataPanel.add(densityDisplayTextBox.graphic(),vertGBC);
        dataPanel.add(temperatureDisplayTextBox.graphic(),vertGBC);
        dataPanel.add(pressureDisplayTextBox.graphic(),vertGBC);

        // Add panels to the control panel
        controlPanel.add(controlsPanel, vertGBC);
        controlPanel.add(setupPanel, vertGBC);
        plotPanel.add(displayCycles.graphic(),vertGBC);
        plotPanel.add(dataPanel, vertGBC);

        //
	    // Configuration tabbed page
        //

	    //slider for scale of display
	    ModifierFunctionWrapper scaleModulator = new ModifierFunctionWrapper(displayBox, "scale");
	    scaleModulator.setFunction(new etomica.util.Function.Linear(0.01, 0.0));
	    scaleSlider = new DeviceSlider(null, scaleModulator);
	    JPanel scaleSliderPanel = new JPanel(new java.awt.GridLayout(0,1));
	    scaleSliderPanel.add(scaleSlider.graphic());	    
	    scaleSlider.getSlider().addChangeListener(new javax.swing.event.ChangeListener() {
	        public void stateChanged(javax.swing.event.ChangeEvent evt) {
	            displayBox.repaint();
	        }
	    });
	    scaleSlider.setMinimum(10);
	    scaleSlider.setMaximum(100);
	    scaleSlider.getSlider().setValue(100);
	    scaleSlider.getSlider().setMajorTickSpacing(10);
	    scaleSlider.getSlider().setMinorTickSpacing(5);
	    scaleSlider.getSlider().setOrientation(1);
	    scaleSlider.getSlider().setLabelTable(scaleSlider.getSlider().createStandardLabels(10));
		
    	displayBoxPanel = new JPanel(new BorderLayout());
    	displayBoxPanel.add(scaleSliderPanel,BorderLayout.EAST);

        //
	    // Plots tabbed page
        //

        plotD = new DisplayPlot();
        plotT = new DisplayPlot();
        plotP = new DisplayPlot();

        JPanel myPlotPanel = new JPanel(new GridLayout(0, 1));
        myPlotPanel.add(plotD.graphic());
        myPlotPanel.add(plotT.graphic());
        myPlotPanel.add(plotP.graphic());
        plotsPane = new JScrollPane(myPlotPanel);

        // Add plots page to tabbed pane
        tabbedPane.add("Plots", plotsPane);
        
        // Default behavior of a SimulationPanel is to
        // show a single graphic.  Switch to show the
        // tabbed page.
        remove(graphicsPanel);
        add(tabbedPane);


		//add meter and display for current kinetic temperature

		thermometer = new MeterTemperature();

        Thread repainter = new Thread() {
            public void run() {
                while (true) {
                    repaint();
                    try{Thread.sleep(repaintSleep);}
                    catch(InterruptedException e){}
                }
            }
        };
        repainter.start();

        if (doRDF) {
            plotRDF = new DisplayPlot();
            JPanel rdfPanel = new JPanel();
            tabbedPane.add("RDF", rdfPanel);
            plotRDF.setDoLegend(false);
            
            rdfPanel.add(plotRDF.graphic());
        }

        boolean pistonHeld = true;

        controlButtons.setSimulation(pc);
        pc.config.setBoundaryPadding(sigma);
        pc.config.initializeCoordinates(pc.box);

        ((SimulationRestart)controlButtons.getReinitButton().getAction()).setConfiguration(pc.config);
        final IIntegrator integrator = pc.integrator;

        ArrayList dataStreamPumps = controlButtons.getDataStreamPumps();
        
        ((ElementSimple)((AtomTypeLeaf)pc.species.getMoleculeType()).getElement()).setMass(mass);
        int D = pc.getSpace().D();


        tUnit = Kelvin.UNIT;

        int thisSleep = -1;
        if (pc.getSpace().D() == 2) {
            dUnit = new UnitRatio(Mole.UNIT, 
                                    new MKS().area());
            Unit[] units = new Unit[] {Bar.UNIT, new PrefixedUnit(Prefix.NANO, Meter.UNIT)};
            double[] exponents = new double[] {1.0, 1.0};
            pUnit = new CompoundUnit(units, exponents);
            thisSleep = integratorSleep;
        }
        else {
            dUnit = new UnitRatio(Mole.UNIT, Liter.UNIT);
            pUnit = Bar.UNIT;
            thisSleep = integratorSleep3D;
        }
        
        densityDisplayTextBox.setLabel("Density ("+dUnit.symbol()+")");
        pressureDisplayTextBox.setLabel("Pressure ("+pUnit.symbol()+")");

        // set up GUI
        pc.ai.setSleepPeriod(thisSleep);
        pc.integrator.removeAllListeners();

        pc.wallPotential.setLongWall(0,true,true);  // left wall
        pc.wallPotential.setLongWall(0,false,true); // right wall
        // skip top wall
        pc.wallPotential.setLongWall(1,false,false);// bottom wall
        pc.wallPotential.setBox(pc.box);  // so it has a boundary
        
        if (displayBox.graphic() != null) {
            displayBoxPanel.remove(displayBox.graphic());
            tabbedPane.remove(displayBoxPanel);
            tabbedPane.remove(blankPanel);
        }
        if (D == 2) {
            tabbedPane.insertTab(displayBox.getLabel(), null, displayBoxPanel, "", 0);
            displayBox.setPixelUnit(new Pixel(400/pc.box.getBoundary().getDimensions().x(1)));
            displayBox.setBox(pc.box);
            displayBox.setAlign(1,DisplayBox.BOTTOM);
            displayBox.canvas.setDrawBoundary(DisplayCanvasInterface.DRAW_BOUNDARY_NONE);
            displayBox.getDrawables().clear();
            displayBox.addDrawable(pc.pistonPotential);
            displayBox.addDrawable(pc.wallPotential);
            displayBoxPanel.add(displayBox.graphic(),BorderLayout.WEST);
            tabbedPane.setSelectedComponent(displayBoxPanel);
        } else {
            tabbedPane.add("Run Faster", blankPanel);
        }
        scaleSlider.setController(pc.controller);

        //  control panel
        ModifierBoolean fixPistonModulator = new ModifierBoolean() {
            public void setBoolean(boolean b) {
                pc.pistonPotential.setStationary(b);
                pressureSlider.getSlider().setEnabled(!b);
                pressureSlider.getTextField().setEnabled(!b);
            }
            public boolean getBoolean() {
                return pc.pistonPotential.isStationary();
            }
        };
        fixPistonButton.setController(pc.controller);
        fixPistonButton.setModifier(fixPistonModulator, "Release piston", "Hold piston");
        fixPistonButton.setPostAction(new ActionPistonUpdate(pc.integrator));
        fixPistonButton.setState(pistonHeld);

        meterCycles = new DataSourceCountTime(pc.integrator);
        displayCycles.setPrecision(6);
        DataPump pump= new DataPump(meterCycles,displayCycles);
        pc.integrator.addIntervalAction(pump);
        displayCycles.setLabel("Simulation time");
        
        //  state panel
        pc.integrator.setIsothermal(tempSlider.isIsothermal());
        pc.integrator.setTemperature(tUnit.toSim(tempSlider.getTemperature()));
        tempSlider.setUnit(tUnit);
        tempSlider.setModifier(new ModifierGeneral(pc.integrator,"temperature"));
        tempSlider.setSliderPostAction(new IntegratorReset(pc.integrator,true));

        potentialSW = new P2SquareWell(pc.getSpace(),sigma,lambda,epsilon,true);
        potentialHS = new P2HardSphere(pc.getSpace(),sigma,true);
        potentialIdeal = new P2Ideal(pc.getSpace());
        
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
        sigModifier.setValue(sigma);
        ModifierGeneral epsModifier = new ModifierGeneral(potentialSW, "epsilon");
        ModifierGeneral lamModifier = new ModifierGeneral(potentialSW, "lambda");
        ModifierGeneral massModifier = new ModifierGeneral(((AtomTypeLeaf)pc.species.getMoleculeType()).getElement(),"mass");
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
        pressureSliderPanel.setBorder(new TitledBorder("Set Pressure ("+pUnit.symbol()+")"));
        Dimension pDim = Pressure.dimension(D);
        pc.pistonPotential.setPressure(pUnit.toSim(pressureSlider.getValue()));
        pressureSlider.setModifier(new ModifierPistonPressure(pc.pistonPotential,pDim));
        pressureSlider.setPostAction(new ActionPistonUpdate(pc.integrator));
        pressureSlider.setController(pc.getController());
        pressureSlider.getSlider().setEnabled(!pc.pistonPotential.isStationary());
        pressureSlider.getTextField().setEnabled(!pc.pistonPotential.isStationary());


        nSlider.setController(pc.getController());
        nSlider.setResetAction(controlButtons.getReinitButton().getAction());
        nSlider.setBox(pc.box);
        nSlider.setSpecies(pc.species);
        nSlider.setMinimum(0);
        nSlider.setMaximum(200);


        if(showTimeControls == true) {
	        doSleepSlider.setModifier(new Modifier() {
	            public String getLabel() {return "";}
	            public Dimension getDimension() {return Undefined.DIMENSION;}
	            public void setValue(double v) {
	                int mySleep = (int)v;
	                if (pc.getSpace().D() == 3) {
	                    integratorSleep3D = mySleep;
	                }
	                else {
	                    integratorSleep = mySleep;
	                }
	                pc.ai.setSleepPeriod(mySleep);
	            }
	            public double getValue() {return integratorSleep;}
	        });

	        pc.integrator.setTimeStep(integratorTimeStepSlider.getValue());
	        integratorTimeStepSlider.setModifier(new Modifier() {
	            public String getLabel() {return "";}
	            public Dimension getDimension() {return Time.DIMENSION;}
	            public void setValue(double v) {
	                pc.integrator.setTimeStep(v);
	            }
	            public double getValue() {return pc.integrator.getTimeStep();}
	        });
	        doSleepSlider.setValue(thisSleep);
        }

        //  data panel
        // the data channel sending the DataTable should be cut off (and hopefully garbage collected)
        plotD.getDataSet().reset();
        plotT.getDataSet().reset();
        plotP.getDataSet().reset();
        
        thermometer.setBox(pc.box);
        AccumulatorHistory temperatureHistory = new AccumulatorHistory();
        temperatureHistory.setTimeDataSource(meterCycles);
        temperatureHistory.getHistory().setHistoryLength(historyLength);
        final AccumulatorAverage temperatureAvg = new AccumulatorAverageCollapsing(100);
        temperatureAvg.setPushInterval(10);
        pump = new DataPump(thermometer,new DataFork(new DataSink[]{temperatureHistory,temperatureAvg}));
        dataStreamPumps.add(pump);
        pc.integrator.addIntervalAction(pump);
        pc.integrator.setActionInterval(pump, dataInterval);
        temperatureHistory.addDataSink(plotT.getDataSet().makeDataSink());
        plotT.setLegend(new DataTag[]{thermometer.getTag()}, "measured");
        temperatureDisplayTextBox.setAccumulator(temperatureAvg);
        temperatureDisplayTextBox.setUnit(tUnit);
        
        DataSource targetTemperatureDataSource = new DataSourceScalar("Target Temperature",Temperature.DIMENSION) {
            public double getDataAsScalar() {
                return tempSlider.getModifier().getValue();
            }
        };
        AccumulatorHistory targetTemperatureHistory = new AccumulatorHistory();
        targetTemperatureHistory.setTimeDataSource(meterCycles);
        targetTemperatureHistory.getHistory().setHistoryLength(historyLength);
        DataPump targetTemperatureDataPump = new DataPump(targetTemperatureDataSource, targetTemperatureHistory);
        pc.integrator.addIntervalAction(targetTemperatureDataPump);
        pc.integrator.setActionInterval(targetTemperatureDataPump, dataInterval);
        targetTemperatureHistory.addDataSink(plotT.getDataSet().makeDataSink());
        plotT.setLegend(new DataTag[]{targetTemperatureDataSource.getTag()}, "target");
        dataStreamPumps.add(targetTemperatureDataPump);

        pressureMeter = new DataSourceWallPressure(pc.getSpace(),pc.pistonPotential);
        pressureMeter.setIntegrator(pc.integrator);
        AccumulatorHistory pressureHistory = new AccumulatorHistory();
        pressureHistory.setTimeDataSource(meterCycles);
        pressureHistory.getHistory().setHistoryLength(historyLength);
        final AccumulatorAverage pressureAvg = new AccumulatorAverageCollapsing(100);
        pressureAvg.setPushInterval(10);
        pump = new DataPump(pressureMeter, new DataFork(new DataSink[]{pressureHistory,pressureAvg}));
        dataStreamPumps.add(pump);
        pc.integrator.addIntervalAction(pump);
        pc.integrator.setActionInterval(pump, dataInterval);
        pressureHistory.addDataSink(plotP.getDataSet().makeDataSink());
        plotP.setLegend(new DataTag[]{pressureMeter.getTag()}, "measured");
        pressureDisplayTextBox.setAccumulator(pressureAvg);
        pressureDisplayTextBox.setUnit(pUnit);

        DataSource targetPressureDataSource = new DataSourceScalar("Target Pressure",Pressure.DIMENSION) {
            public double getDataAsScalar() {
                return pUnit.toSim(pressureSlider.getValue());
            }
        };
        AccumulatorHistory targetPressureHistory = new AccumulatorHistory();
        targetPressureHistory.setTimeDataSource(meterCycles);
        targetPressureHistory.getHistory().setHistoryLength(historyLength);
        pump = new DataPump(targetPressureDataSource, targetPressureHistory);
        pc.integrator.addIntervalAction(pump);
        pc.integrator.setActionInterval(pump, dataInterval);
        targetPressureHistory.addDataSink(plotP.getDataSet().makeDataSink());
        plotP.setLegend(new DataTag[]{targetPressureDataSource.getTag()}, "target");
        dataStreamPumps.add(pump);

        densityMeter = new MeterDensity(pc.getSpace()); //pc.pistonPotential,1);
        densityMeter.setBox(pc.box);
        AccumulatorHistory densityHistory = new AccumulatorHistory();
        densityHistory.setTimeDataSource(meterCycles);
        densityHistory.getHistory().setHistoryLength(historyLength);
        final AccumulatorAverage densityAvg = new AccumulatorAverageCollapsing(100);
        densityAvg.setPushInterval(10);
        pump = new DataPump(densityMeter,new DataFork(new DataSink[]{densityAvg, densityHistory}));
        dataStreamPumps.add(pump);
        pc.integrator.addIntervalAction(pump);
        pc.integrator.setActionInterval(pump, dataInterval);
        densityHistory.addDataSink(plotD.getDataSet().makeDataSink());
        plotD.setLegend(new DataTag[]{densityMeter.getTag()}, "measured");
        densityDisplayTextBox.setAccumulator(densityAvg);
        densityDisplayTextBox.setUnit(dUnit);
        
        plotD.setUnit(dUnit);
        plotT.setUnit(tUnit);
        plotP.setUnit(pUnit);

        plotP.getPlot().setTitle("Pressure ("+pUnit.symbol()+")");
        plotT.getPlot().setTitle("Temperature ("+tUnit.symbol()+")");
        plotD.getPlot().setTitle("Density ("+dUnit.symbol()+")");
 
        // Set the size of the plots and the scoll pane containing the plots.
        // Want 2 of the 3 plots displayed
        java.awt.Dimension d = plotT.getPlot().getPreferredSize();
        d.height = 230;
        plotT.getPlot().setSize(d);
        plotP.getPlot().setSize(d);
        plotD.getPlot().setSize(d);

        d.width += 40;
        d.height = d.height * 2 + 40;
        plotsPane.setPreferredSize(d);

        if (doRDF) {
            plotRDF.getDataSet().reset();
            double rdfCutoff = 10;
            final MeterRDFCylinder meterRDF = new MeterRDFCylinder(pc.getSpace());
            meterRDF.setBox(pc.box);
            meterRDF.getXDataSource().setXMax(rdfCutoff);
            meterRDF.setPotential(pc.pistonPotential);
            pump = new DataPump(meterRDF, plotRDF.getDataSet().makeDataSink());
            pc.integrator.addIntervalAction(meterRDF);
            dataStreamPumps.add(pump);
            pc.integrator.addIntervalAction(pump);
            
            controlButtons.getResetAveragesButton().setPostAction(new Action() {
                public void actionPerformed() {
                    meterRDF.reset();
                }
            });
        }
        
        if (doConfigButton) {
            configButton.setController(pc.getController());
            configButton.setAction(new ActionConfigWindow(pc.box));

            velocityButton.setController(pc.getController());
            velocityButton.setAction(new ActionVelocityWindow(pc.box));
        }

        // re-initialize the integrator expclitly.  This resets the piston
        // position back to the top.  Also, force the data into the
        // display boxes and repaint the boxes.
        controlButtons.getReinitButton().setPostAction(new Action() {
            public void actionPerformed() {
                try {
                    integrator.initialize();

                    densityDisplayTextBox.putData(densityAvg.getData());
                    densityDisplayTextBox.repaint();

                    temperatureDisplayTextBox.putData(temperatureAvg.getData());
                    temperatureDisplayTextBox.repaint();

                    pressureDisplayTextBox.putData(pressureAvg.getData());
                    pressureDisplayTextBox.repaint();
                }
                catch (ConfigurationOverlapException e) {
                    throw new RuntimeException(e);
                }
            }
        });

    }
    
    public void setPotential(String potentialDesc) {
        final boolean HS = potentialDesc.equals("Repulsion only"); 
        final boolean SW = potentialDesc.equals("Repulsion and attraction"); 
        pc.controller.doActionNow( new Action() {
            public void actionPerformed() {
                if (HS) {
                    potentialHS.setBox(pc.box);
                    pc.potentialWrapper.setPotential(potentialHS);
                }
                else if (SW) {
                    potentialSW.setBox(pc.box);
                    pc.potentialWrapper.setPotential(potentialSW);
                }
                else {
                    potentialIdeal.setBox(pc.box);
                    pc.potentialWrapper.setPotential(potentialIdeal);
                }
                try {
                    if (pc.integrator.isInitialized()) {
                        pc.integrator.reset();
                    }
                } catch(ConfigurationOverlapException e) {}
            }
        });
    }

    protected class ModifierAtomDiameter implements Modifier {

        public void setValue(double d) {
            //assume one type of atom
            ((AtomTypeSphere)pc.species.getMoleculeType()).setDiameter(d);
            PistonCylinderGraphic.this.potentialHS.setCollisionDiameter(d);
            PistonCylinderGraphic.this.potentialSW.setCoreDiameter(d);
            pc.pistonPotential.setCollisionRadius(0.5*d);
            pc.wallPotential.setCollisionRadius(0.5*d);
            sigma = d;
            displayBox.repaint();
            pc.config.setBoundaryPadding(sigma);
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

    public static void main(String[] args) {
        PistonCylinderGraphic pcg = new PistonCylinderGraphic();
        pcg.init();
		SimulationGraphic.makeAndDisplayFrame(pcg, APP_NAME);
    }

    public static class Applet extends javax.swing.JApplet {

        public void init() {
            PistonCylinderGraphic pcg = new PistonCylinderGraphic();
            String doConfigButtonStr = getParameter("doConfigButton");
            if (doConfigButtonStr != null) {
                pcg.setDoConfigButton(Boolean.valueOf(doConfigButtonStr).booleanValue());
            }
            String doRDFStr = getParameter("doRDF");
            if (doRDFStr != null) {
                pcg.setDoRDF(Boolean.valueOf(doRDFStr).booleanValue());
            }
            String doNSelectorStr = getParameter("doRDF");

            pcg.init();
            getContentPane().add(pcg);
        }

        private static final long serialVersionUID = 1L;
    }
}
