package etomica.modules.pistoncylinder;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.event.ItemListener;
import java.util.ArrayList;

import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.border.TitledBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import etomica.action.Action;
import etomica.action.ActionGroupSeries;
import etomica.action.IntegratorReset;
import etomica.action.SimulationRestart;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IVector;
import etomica.atom.AtomSet;
import etomica.atom.AtomTypeSphere;
import etomica.atom.IAtomPositioned;
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
import etomica.graphics.DeviceDelaySlider;
import etomica.graphics.DeviceNSelector;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DeviceThermoSlider;
import etomica.graphics.DeviceToggleButton;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.DisplayCanvasInterface;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.graphics.DeviceBox.LabelType;
import etomica.integrator.IIntegrator;
import etomica.modifier.Modifier;
import etomica.modifier.ModifierBoolean;
import etomica.modifier.ModifierFunctionWrapper;
import etomica.modifier.ModifierGeneral;
import etomica.potential.P2HardSphere;
import etomica.potential.P2Ideal;
import etomica.potential.P2SquareWell;
import etomica.potential.Potential2HardSphericalWrapper;
import etomica.space.Space;
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
import etomica.units.Unit;
import etomica.units.UnitRatio;
import etomica.units.systems.MKS;


public class PistonCylinderGraphic extends SimulationGraphic {
    
	private static final String APP_NAME = "Piston Cylinder";
	private final static int REPAINT_INTERVAL = 1;

    public PistonCylinder pc;
    public Potential2HardSphericalWrapper potentialWrapper;
    public P2HardSphere potentialHS;
    public P2SquareWell potentialSW;
    public P2Ideal potentialIdeal;
    public DataSourceCountTime meterCycles;
    public DisplayTextBox displayCycles;
    public MeterTemperature thermometer;
    public DisplayBox displayBox;
    public DeviceButton configButton, velocityButton;
    public DeviceButton goFastButton;
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
    public DeviceBox densityBox;
	protected DisplayTextBoxesCAE densityDisplayTextBox, temperatureDisplayTextBox, pressureDisplayTextBox;
    public JPanel blankPanel = new JPanel();
    public JScrollPane plotsPane;
    public int historyLength;
    public DataSourceWallPressure pressureMeter;
    public int dataInterval;
    public Unit eUnit;
    public double lambda, epsilon, mass, sigma;
    public DeviceSlider doSleepSlider, integratorTimeStepSlider;

    protected boolean doConfigButton = false;
    protected boolean doRDF = false;
    protected boolean doDensityInput = false;
    protected boolean doNMoleculeSlider = false;

    /**
     * Creates a PistonCylinder graphic instance.  init() must be called before
     * this can be used.
     */
    public PistonCylinderGraphic(PistonCylinder sim, Space _space) {
    	super(sim, TABBED_PANE, APP_NAME, REPAINT_INTERVAL, _space);
    	pc = sim;
    }

    /**
     * Enable the config (and velocity) buttons.  This must be called before
     * init() is called.
     */
    public void setDoConfigButton(boolean newDoConfigButton) {
        doConfigButton = newDoConfigButton;
    }
    
    /**
     * Enable the RDF plot.  This must be called before init() is called.
     */
    public void setDoRDF(boolean newDoRDF) {
        doRDF = newDoRDF;
    }
    
    /**
     * Enable the RDF plot.  This must be called before init() is called.
     */
    public void setDoDensityInput(boolean newDoDensityInput) {
        doDensityInput = newDoDensityInput;
    }

    /**
     * Enable the # of molecules slider.  This must be called before init() is called.
     */
    public void setDoNMoleculeSlider(boolean newDoNMoleculeSlider) {
        doNMoleculeSlider = newDoNMoleculeSlider;
    }

    public void setRepaintInterval(int newRepaintInterval) {
        pc.integrator.setActionInterval(getPaintAction(pc.box), newRepaintInterval);
    }
    
    /**
     * Initialize all the bits based on previously set doConfigButton and doRDF.
     */
    public void init() {

        displayBox = getDisplayBox(pc.box);
        displayBox.setColorScheme(new ColorSchemeByType());
        if (pc.getSpace().D() == 3) {
            pc.integrator.setActionInterval(getPaintAction(pc.box), 1);
            ((DisplayBoxCanvasG3DSys)displayBox.canvas).addPlane(new PistonPlane(pc.pistonPotential));
        }

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

        // Release/Hold piston button
        fixPistonButton = new DeviceToggleButton(null);
        ((JPanel)getController().graphic()).add(fixPistonButton.graphic());

        // Simulation Time
        displayCycles = new DisplayTextBox();

        // Show buttons (added to panel 
        if (doConfigButton) {
            JPanel configPanel = new JPanel(new GridBagLayout());


            configButton = new DeviceButton(null);
            configButton.setLabel("Show Config");
            configPanel.add(configButton.graphic(),horizGBC);

            velocityButton = new DeviceButton(null);
            velocityButton.setLabel("Show Velocities");
            configPanel.add(velocityButton.graphic(),horizGBC);
            getPanel().controlPanel.add(configPanel,vertGBC);
        }
        
        if (pc.getSpace().D() == 3) {
            goFastButton = new DeviceButton(pc.getController());
            goFastButton.setLabel("Go Fast");
            getPanel().controlPanel.add(goFastButton.graphic(), vertGBC);
        }

        //
        // State tabbed pane page
        //

        tempSlider = new DeviceThermoSlider(pc.getController());
        tempSlider.setShowValues(true);
        tempSlider.setEditValues(true);
        tempSlider.setMinimum(0);
        tempSlider.setMaximum(1000);
        tempSlider.setAdiabatic();
        tempSlider.setSliderMajorValues(4);
        tempSlider.setIntegrator(pc.integrator);
        tempSlider.setTemperature(300);

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

        if (doDensityInput) {
            densityBox = new DeviceBox();
            densityBox.setController(pc.getController());
        }

        if (doNMoleculeSlider) {
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
                    
                    displayBox.repaint();
                }
            });

            nSliderPanel = new JPanel(new GridLayout(0,1));
            nSliderPanel.setBorder(new TitledBorder(null, "Number of Molecules", TitledBorder.CENTER, TitledBorder.TOP));
            nSlider.setShowBorder(false);
            nSlider.setNMajor(4);
            nSliderPanel.add(nSlider.graphic());
        }


        // Add all state page sub panels onto a single panel
        JPanel statePanel = new JPanel(new GridBagLayout());
        GridBagConstraints gbc2 = new GridBagConstraints();
        gbc2.gridx = 0;  gbc2.gridy = 0;
        statePanel.add(tempSlider.graphic(), gbc2);
        gbc2.gridx = 0;  gbc2.gridy = 1;
        statePanel.add(pressureSliderPanel, gbc2);

        if (doNMoleculeSlider) {
            gbc2.gridx = 0;  gbc2.gridy = 2;
            statePanel.add(nSliderPanel, gbc2);
        }

        if (doDensityInput) {
            gbc2.gridx = 0;  gbc2.gridy = 3;
            gbc2.fill = GridBagConstraints.HORIZONTAL;
            statePanel.add(densityBox.graphic(), gbc2);
        }

        //
        // Potential tabbed pane page
        //
        
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


        //
	    // Configuration tabbed page
        //

	    //slider for scale of display
        JPanel scaleSliderPanel = null;
        if (pc.getSpace().D() == 2) {
            scaleSliderPanel = new JPanel();
    	    ModifierFunctionWrapper scaleModulator = new ModifierFunctionWrapper(displayBox, "scale");
    	    scaleModulator.setFunction(new etomica.util.Function.Linear(0.01, 0.0));
    	    scaleSlider = new DeviceSlider(null, scaleModulator);
    	    scaleSlider.setShowValues(false);
    	    scaleSliderPanel.setBorder(new TitledBorder(null, "Graphic Size", TitledBorder.CENTER, TitledBorder.TOP));
    	    scaleSliderPanel.add(scaleSlider.graphic());
    	    scaleSlider.getSlider().addChangeListener(new javax.swing.event.ChangeListener() {
    	        public void stateChanged(javax.swing.event.ChangeEvent evt) {
    	        	displayBox.repaint();
    	        }
    	    });
    	    scaleSlider.setMinimum(10);
    	    scaleSlider.setMaximum(100);
    	    scaleSlider.getSlider().setValue(100);
    	    scaleSlider.setNMajor(0);
    	    scaleSlider.setSliderVerticalOrientation(false);
    	    java.util.Hashtable scaleLabels = new java.util.Hashtable();
    	    scaleLabels.put(new Integer(10), new JLabel( "min", JLabel.CENTER ));
    	    scaleLabels.put(new Integer(100), new JLabel( "max", JLabel.CENTER ));
    	    scaleSlider.getSlider().setLabelTable(scaleLabels);
        }

        DeviceDelaySlider delaySlider = new DeviceDelaySlider(pc.getController(), (ActivityIntegrate)pc.getController().getAllActions()[0]);
        
        // Add panels to the control panel
        getPanel().controlPanel.add(setupPanel, vertGBC);
        if (scaleSliderPanel != null) {
            getPanel().controlPanel.add(scaleSliderPanel, vertGBC);
        }
        getPanel().controlPanel.add(delaySlider.graphic(), vertGBC);
        add(displayCycles);
        add(densityDisplayTextBox);
        add(temperatureDisplayTextBox);
        add(pressureDisplayTextBox);

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
        getPanel().tabbedPane.add("Plots", plotsPane);

		//add meter and display for current kinetic temperature

		thermometer = new MeterTemperature(space.D());

        if (doRDF) {
            plotRDF = new DisplayPlot();
            JPanel rdfPanel = new JPanel();
            getPanel().tabbedPane.add("RDF", rdfPanel);
            plotRDF.setDoLegend(false);
            
            rdfPanel.add(plotRDF.graphic());
        }

        boolean pistonHeld = true;

        pc.config.setBoundaryPadding(sigma);
        pc.config.initializeCoordinates(pc.box);

        ((SimulationRestart)getController().getReinitButton().getAction()).setConfiguration(pc.config);
        final IIntegrator integrator = pc.integrator;

        ArrayList dataStreamPumps = getController().getDataStreamPumps();
        
        ((ElementSimple)pc.species.getLeafType().getElement()).setMass(mass);
        int D = pc.getSpace().D();

        tUnit = Kelvin.UNIT;

        if (pc.getSpace().D() == 2) {
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
        
        densityDisplayTextBox.setLabel("Density ("+dUnit.symbol()+")");
        pressureDisplayTextBox.setLabel("Pressure ("+pUnit.symbol()+")");

        pc.wallPotential.setLongWall(0,true,true);  // left wall
        pc.wallPotential.setLongWall(0,false,true); // right wall
        // skip top wall
        pc.wallPotential.setLongWall(1,false,false);// bottom wall
        pc.wallPotential.setBox(pc.box);  // so it has a boundary

/* WILL BE NEEDED IF DIMENSION RADIO BUTTONS USED
        if (displayBox.graphic() != null) {
            remove(displayBox);
            getPanel().tabbedPane.remove(blankPanel);
        }
*/

        if (D == 2) {
            displayBox.setPixelUnit(new Pixel(400/pc.box.getBoundary().getDimensions().x(1)));
        }
        else {
            displayBox.setPixelUnit(new Pixel(40/pc.box.getBoundary().getDimensions().x(1)));
        }
        displayBox.setAlign(1,DisplayBox.BOTTOM);
        displayBox.canvas.setDrawBoundary(DisplayCanvasInterface.DRAW_BOUNDARY_NONE);
        displayBox.getDrawables().clear();
        if (pc.getSpace().D() == 2) {
            // doesn't actually work for 3D
            displayBox.addDrawable(pc.pistonPotential);
            displayBox.addDrawable(pc.wallPotential);
        }
        if (scaleSlider != null) {
            // doesn't actually have any effect for 3D
            scaleSlider.setController(pc.getController());
        }

        //  control panel
        ModifierBoolean fixPistonModulator = new ModifierBoolean() {
            public void setBoolean(boolean b) {
                pc.pistonPotential.setStationary(b);
                pressureSlider.getSlider().setEnabled(!b);
                pressureSlider.getTextField().setEnabled(!b);
                if (doDensityInput) {
                    densityBox.setEditable(b);
                }
            }
            public boolean getBoolean() {
                return pc.pistonPotential.isStationary();
            }
        };
        fixPistonButton.setController(pc.getController());
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
        ModifierGeneral massModifier = new ModifierGeneral(pc.species.getLeafType().getElement(),"mass");
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
        pressureSliderPanel.setBorder(new TitledBorder(null, "Set Pressure ("+pUnit.symbol()+")", TitledBorder.CENTER, TitledBorder.TOP));
        Dimension pDim = Pressure.dimension(D);
        double p = pUnit.toSim(pressureSlider.getValue());
        pc.pistonPotential.setPressure(D == 3 ? -p : p);
        pressureSlider.setModifier(new ModifierPistonPressure(pc.pistonPotential,pDim));
        pressureSlider.setPostAction(new ActionPistonUpdate(pc.integrator));
        pressureSlider.setController(pc.getController());
        pressureSlider.getSlider().setEnabled(!pc.pistonPotential.isStationary());
        pressureSlider.getTextField().setEnabled(!pc.pistonPotential.isStationary());


        if (doNMoleculeSlider) {
            nSlider.setController(pc.getController());
            nSlider.setResetAction(new Action() {
                public void actionPerformed() {
                    pc.integrator.resetPiston();
                    getController().getReinitButton().getAction().actionPerformed();
                }
            });
            nSlider.setBox(pc.box);
            nSlider.setSpecies(pc.species);
            nSlider.setMinimum(0);
            nSlider.setMaximum(200);
        }

        if (doDensityInput) {
            densityBox.setLabelType(LabelType.BORDER);
            densityBox.setLabel("Set Density ("+dUnit.symbol()+")");
            densityBox.setUnit(dUnit);
            densityBox.setModifier(new ModifierPistonDensity());
            densityBox.setPostAction(new ActionGroupSeries(new Action[]{new IntegratorReset(pc.integrator,true), getPaintAction(pc.box)}));
        }

        //  data panel
        // the data channel sending the DataTable should be cut off (and hopefully garbage collected)
        plotD.getDataSet().reset();
        plotT.getDataSet().reset();
        plotP.getDataSet().reset();
        
        thermometer.setBox(pc.box);
        final AccumulatorHistory temperatureHistory = new AccumulatorHistory();
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
        final AccumulatorHistory targetTemperatureHistory = new AccumulatorHistory();
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
        final AccumulatorHistory pressureHistory = new AccumulatorHistory();
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
        final AccumulatorHistory targetPressureHistory = new AccumulatorHistory();
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
        final AccumulatorHistory densityHistory = new AccumulatorHistory();
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
            
            getController().getResetAveragesButton().setPostAction(new Action() {
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
        getController().getReinitButton().setPostAction(new Action() {
            public void actionPerformed() {
                pc.integrator.resetPiston();
                if (pc.integrator.isInitialized()) {
                    try {
                        pc.integrator.reset();
                    }
                    catch (ConfigurationOverlapException e) {}
                }

                densityDisplayTextBox.putData(densityAvg.getData());
                densityDisplayTextBox.repaint();

                temperatureDisplayTextBox.putData(temperatureAvg.getData());
                temperatureDisplayTextBox.repaint();

                pressureDisplayTextBox.putData(pressureAvg.getData());
                pressureDisplayTextBox.repaint();
                
                displayBox.repaint();
                
                displayCycles.putData(meterCycles.getData());
                displayCycles.repaint();
            }
        });
        
        if (pc.getSpace().D() == 3) {
            goFastButton.setAction(new Action() {
                public void actionPerformed() {
                    if (isFast) {
                        isFast = false;
                        goFastButton.setLabel("Go Fast");
                        pc.integrator.setTimeStep(1);
                        setRepaintInterval(1);
                        pc.integrator.setEventInterval(1);
                        densityHistory.setActive(true);
                        temperatureHistory.setActive(true);
                        pressureHistory.setActive(true);
                        targetTemperatureHistory.setActive(true);
                        targetPressureHistory.setActive(true);
                    }
                    else {
                        isFast = true;
                        goFastButton.setLabel("Go Slower");
                        pc.integrator.setTimeStep(10);
                        pc.ai.setSleepPeriod(0);
                        setRepaintInterval(10000);
                        pc.integrator.setEventInterval(1);
                        densityHistory.setActive(false);
                        temperatureHistory.setActive(false);
                        pressureHistory.setActive(false);
                        targetTemperatureHistory.setActive(false);
                        targetPressureHistory.setActive(false);
                    }
                }
                
                protected boolean isFast = false;
            });
        }
    }
    
    public void setPotential(String potentialDesc) {
        final boolean HS = potentialDesc.equals("Repulsion only"); 
        final boolean SW = potentialDesc.equals("Repulsion and attraction"); 
        pc.getController().doActionNow( new Action() {
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

    public class ModifierPistonDensity implements Modifier {
        public Dimension getDimension() {
            return dUnit.dimension();
        }

        public String getLabel() {
            return "density";
        }

        public double getValue() {
            return pc.box.getDensity();
        }

        public void setValue(double newValue) {
            double oldDensity = pc.box.getDensity();
            double maxDensity;
            int D = pc.getSpace().D();
            if (D == 2) {
                maxDensity = 2.0 / (Math.sqrt(3) * sigma * sigma);
            }
            else {
                maxDensity = Math.sqrt(2) / (sigma*sigma*sigma);
            }
            if (newValue > maxDensity && newValue > oldDensity) {
                if (oldDensity > maxDensity) {
                    return;
                }
                newValue = maxDensity;
            }
            IVector boxDim = pc.box.getBoundary().getDimensions();
            AtomSet leafList = pc.box.getLeafList();
            double yShift = 0.5*(boxDim.x(1)-sigma);
            if (D == 2) {
                yShift = -yShift;
            }
            if (newValue > oldDensity) {
                // scale atom positions
                for (int i=0; i<leafList.getAtomCount(); i++) {
                    IVector pos = ((IAtomPositioned)leafList.getAtom(i)).getPosition();
                    double y = (pos.x(1)+yShift) * (oldDensity / newValue) - yShift;
                    pos.setX(1, y);
                }
            }
            yShift += (D==2 ? 1:-1) * 0.5*sigma;
            double pistonY = pc.pistonPotential.getWallPosition();
            pistonY = (pistonY + yShift) * (oldDensity / newValue) - yShift;
            pc.pistonPotential.setWallPosition(pistonY);
        }
    }

    protected class ModifierAtomDiameter implements Modifier {

        public void setValue(double d) {
            //assume one type of atom
            ((AtomTypeSphere)pc.species.getLeafType()).setDiameter(d);
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
        PistonCylinder sim = new PistonCylinder(2);
        PistonCylinderGraphic pcg = new PistonCylinderGraphic(sim, sim.getSpace());
        pcg.setDoRDF(true);
        pcg.setDoDensityInput(true);
        pcg.setDoConfigButton(true);
        pcg.setDoNMoleculeSlider(true);
        pcg.init();
		SimulationGraphic.makeAndDisplayFrame(pcg.getPanel(), APP_NAME);
    }

    public static class Applet extends javax.swing.JApplet {

        public void init() {
            int D = 2;
            String dimStr = getParameter("dim");
            if (dimStr != null) {
                D = Integer.parseInt(dimStr);
            }
            PistonCylinder sim = new PistonCylinder(D);
            PistonCylinderGraphic pcg = new PistonCylinderGraphic(sim, sim.getSpace());
            String doConfigButtonStr = getParameter("doConfigButton");
            if (doConfigButtonStr != null) {
                pcg.setDoConfigButton(Boolean.valueOf(doConfigButtonStr).booleanValue());
            }
            String doRDFStr = getParameter("doRDF");
            if (doRDFStr != null) {
                pcg.setDoRDF(Boolean.valueOf(doRDFStr).booleanValue());
            }
            String doDensityInputStr = getParameter("doDensityInput");
            if (doDensityInputStr != null) {
                pcg.setDoDensityInput(Boolean.valueOf(doDensityInputStr).booleanValue());
            }
            String doNMoleculeSlider = getParameter("doNMoleculeSlider");
            if (doNMoleculeSlider != null) {
                pcg.setDoNMoleculeSlider(Boolean.valueOf(doNMoleculeSlider).booleanValue());
            }
            pcg.init();
            String repaintIntervalStr = getParameter("repaintInterval");
            if (repaintIntervalStr != null) {
                pcg.setRepaintInterval(Integer.valueOf(repaintIntervalStr).intValue());
            }
            getContentPane().add(pcg.getPanel());
        }

        private static final long serialVersionUID = 1L;
    }
}
