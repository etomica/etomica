/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.osmosis;

import etomica.action.IAction;
import etomica.action.SimulationRestart;
import etomica.action.activity.ActivityIntegrate;
import etomica.config.ConfigurationLatticeWithPlane;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.DataPump;
import etomica.data.DataSourceCountTime;
import etomica.data.meter.MeterLocalMoleFraction;
import etomica.data.meter.MeterTemperature;
import etomica.graphics.*;
import etomica.integrator.IntegratorListenerAction;
import etomica.lattice.LatticeCubicSimple;
import etomica.math.geometry.Cuboid;
import etomica.math.geometry.Plane;
import etomica.math.geometry.Rectangle;
import etomica.potential.P1HardBoundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space2d.Space2D;
import etomica.space2d.Vector2D;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.units.Kelvin;
import etomica.units.Unit;
import etomica.util.Constants.CompassDirection;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.awt.*;
import java.util.ArrayList;

/**
 * Osmosis module.
 * @author Jhumpa Adhikari
 * @author Andrew Schultz
 */

public class Osmosis extends SimulationGraphic {

	private final static String APP_NAME = "Osmosis";
	private final static int REPAINT_INTERVAL = 1;

    public DataSourceCountTime cycles;
    public DisplayTextBox displayCycles;
    public MeterOsmoticPressure osmosisPMeter;
    public MeterLocalMoleFraction moleFractionRight, moleFractionLeft;
    public OsmosisSim sim;
    private DeviceThermoSlider temperatureSelect;

    public Osmosis(OsmosisSim simulation, Space _space) {

    	super(simulation, TABBED_PANE, APP_NAME, REPAINT_INTERVAL);

        ArrayList<DataPump> dataStreamPumps = getController().getDataStreamPumps();
        
    	sim = simulation;
        final int thickness = 4;

        Unit tUnit = Kelvin.UNIT;

        ConfigurationLatticeWithPlane config = null;

	    //display of box
        final DisplayBox displayBox = getDisplayBox(sim.box);
        ColorSchemeByType colorScheme = new ColorSchemeByType();

        colorScheme.setColor(sim.speciesSolvent.getLeafType(), Color.blue);
        colorScheme.setColor(sim.speciesSolute.getLeafType(), Color.red);
        displayBox.setColorScheme(colorScheme);
        displayBox.setAlign(1,DisplayBox.CENTER);
        displayBox.setOriginShift(0, thickness);
        displayBox.setOriginShift(1, -thickness);
        if (sim.getSpace() instanceof Space2D) {
            displayBox.addDrawable(new MyWall());
            config = new ConfigurationLatticeWithPlane(new LatticeCubicSimple(space, 1.0), null, space);
        }
        else if (sim.getSpace() instanceof Space3D) {
        	Plane plane = new Plane(sim.getSpace());
        	((etomica.graphics.DisplayBoxCanvasG3DSys)displayBox.canvas).addPlane(plane);
            config = new ConfigurationLatticeWithPlane(new LatticeCubicSimple(space, 1.0), plane, space); 
            config.addSpecies(sim.speciesSolvent);
            config.addSpecies(sim.speciesSolute);
        }

        config.initializeCoordinates(sim.box);

        SimulationRestart simRestart = getController().getSimRestart();
        final IAction resetAction = simRestart.getDataResetAction();
        
        simRestart.setConfiguration(config);

        cycles = new DataSourceCountTime();
        displayCycles = new DisplayTimer(sim.integrator);
        displayCycles.setLabel("Cycle Time");
	    displayCycles.setPrecision(6);	

	    // Right side of membrane osmotic
        osmosisPMeter = new MeterOsmoticPressure(new P1HardBoundary[]{sim.boundaryHardA,sim.boundaryHardB}, sim.integrator);
        final AccumulatorAverageCollapsing osmosisPMeterAvg = new AccumulatorAverageCollapsing();
        final DataPump osmosisPump = new DataPump(osmosisPMeter, osmosisPMeterAvg);
        dataStreamPumps.add(osmosisPump);
        IntegratorListenerAction osmosisPumpListener = new IntegratorListenerAction(osmosisPump);
        sim.integrator.getEventManager().addListener(osmosisPumpListener);
        osmosisPumpListener.setInterval(40);
        sim.integrator.setTimeStep(0.01);
        final DisplayTextBoxesCAE osmoticBox = new DisplayTextBoxesCAE();
        osmoticBox.setAccumulator(osmosisPMeterAvg);
        osmoticBox.setPrecision(5);
//        osmoticBox.setLabelType(LabelType.BORDER);
        osmoticBox.setLabel("Osmotic Pressure (PV/Nk)");

        //
        // Delay panel
        //
        DeviceDelaySlider delaySlider = new DeviceDelaySlider(sim.getController());
        
        //
        // temperature panel
        //

	    temperatureSelect = new DeviceThermoSlider(sim.getController(), sim.integrator);
	    temperatureSelect.setUnit(tUnit);
	    temperatureSelect.setMaximum(1000);
	    temperatureSelect.setTemperature(300);
	    temperatureSelect.setIsothermal();
	    temperatureSelect.setSliderPostAction(resetAction);
	    temperatureSelect.setRadioGroupPostAction(resetAction);
		MeterTemperature thermometer = new MeterTemperature(sim.box, space.D());
		DisplayTextBox tBox = new DisplayTextBox();
        DataPump tempPump = new DataPump(thermometer, tBox);
        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(tempPump));
		tBox.setUnit(tUnit);
		tBox.setLabel("Measured Temperature");
		tBox.setLabelPosition(CompassDirection.NORTH);
		tBox.setPrecision(3);

        // Right side of membrane mole fraction
        moleFractionRight = new MeterLocalMoleFraction(space, sim.box);
        Vector dimensions = sim.box.getBoundary().getBoxSize();

        if (sim.getSpace() instanceof Space2D) { // 2D
            moleFractionRight.setShape(new Rectangle(sim.getSpace(), dimensions.getX(0)*0.5, dimensions.getX(1)));
            moleFractionRight.setShapeOrigin(new Vector2D(dimensions.getX(0)*0.25, 0));
        }
        else if (sim.getSpace() instanceof Space3D) { // 3D
            moleFractionRight.setShape(new Cuboid(sim.getSpace(), dimensions.getX(0)*0.5, dimensions.getX(1), dimensions.getX(2)));
            moleFractionRight.setShapeOrigin(new Vector3D(dimensions.getX(0)*0.25, 0, 0));
        }

        moleFractionRight.setSpecies(sim.speciesSolute);
        final AccumulatorAverageCollapsing moleFractionAvgRight = new AccumulatorAverageCollapsing();
        final DataPump molePumpRight = new DataPump(moleFractionRight, moleFractionAvgRight);
        dataStreamPumps.add(molePumpRight);
        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(molePumpRight));
        final DisplayTextBoxesCAE rightMFBox = new DisplayTextBoxesCAE();
        rightMFBox.setAccumulator(moleFractionAvgRight);
        rightMFBox.setPrecision(5);
//        rightMFBox.setLabelType(LabelType.BORDER);
        rightMFBox.setLabel("Mole Fraction Right of Membrane(nSolute/nSolution)");

        // Left side of membrane mole fraction
        moleFractionLeft = new MeterLocalMoleFraction(space, sim.box);

        if (sim.getSpace() instanceof Space2D) { // 2D
            moleFractionLeft.setShape(new Rectangle(sim.getSpace(), dimensions.getX(0)*0.5, dimensions.getX(1)));
            moleFractionLeft.setShapeOrigin(new Vector2D(-dimensions.getX(0)*0.25, 0));
        }
        else if (sim.getSpace() instanceof Space3D) { // 3D
            moleFractionLeft.setShape(new Cuboid(sim.getSpace(), dimensions.getX(0)*0.5, dimensions.getX(1), dimensions.getX(2)));

            moleFractionLeft.setShapeOrigin(new Vector3D(-dimensions.getX(0)*0.25, 0, 0));
        }

        moleFractionLeft.setSpecies(sim.speciesSolute);
        final AccumulatorAverageCollapsing moleFractionAvgLeft = new AccumulatorAverageCollapsing();
        final DataPump molePumpLeft = new DataPump(moleFractionLeft, moleFractionAvgLeft);
        dataStreamPumps.add(molePumpLeft);
        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(molePumpLeft));
        final DisplayTextBoxesCAE leftMFBox = new DisplayTextBoxesCAE();
        leftMFBox.setAccumulator(moleFractionAvgLeft);
        leftMFBox.setPrecision(5);
//        leftMFBox.setLabelType(LabelType.BORDER);
        leftMFBox.setLabel("Mole Fraction Left of Membrane(nSolute/nSolution)");

        ChangeListener cl = new ChangeListener() {
        	public void stateChanged(ChangeEvent evt) {
				molePumpLeft.actionPerformed();
				leftMFBox.putData(moleFractionAvgLeft.getData());

				molePumpRight.actionPerformed();
				rightMFBox.putData(moleFractionAvgRight.getData());
				osmoticBox.putData(osmosisPMeterAvg.getData());
				getDisplayBox(sim.box).graphic().repaint();
        	}
        };

        InitializeMolecules initPanel = new InitializeMolecules(config, simRestart);
        initPanel.addStateChangedListener(cl);

        // left side panel for mole fraction
        JPanel leftMoleFractionPanel = new JPanel(new FlowLayout());
		leftMoleFractionPanel.setBorder(new TitledBorder(null, "Mole Fraction (nSolute/nSolution)", TitledBorder.CENTER, TitledBorder.TOP));
		leftMoleFractionPanel.add(leftMFBox.graphic());

        // left side metrics panel
        JPanel leftMetricsPanel = new JPanel(new GridLayout(0, 1));
        leftMetricsPanel.setBorder(new TitledBorder(null, "Left of Membrane", TitledBorder.CENTER, TitledBorder.TOP));
        leftMetricsPanel.add(leftMoleFractionPanel);

        // right side panel for mole fraction
        JPanel rightMoleFractionPanel = new JPanel(new FlowLayout());
		rightMoleFractionPanel.setBorder(new TitledBorder(null, "Mole Fraction (nSolute/nSolution)", TitledBorder.CENTER, TitledBorder.TOP));
		rightMoleFractionPanel.add(rightMFBox.graphic());

        // right side metrics panel
        JPanel rightMetricsPanel = new JPanel(new GridLayout(0, 1));
        rightMetricsPanel.setBorder(new TitledBorder(null, "Right of Membrane", TitledBorder.CENTER, TitledBorder.TOP));
        rightMetricsPanel.add(rightMoleFractionPanel);

        GridBagConstraints vertGBC = SimulationPanel.getVertGBC();

        add(displayCycles);
        add(tBox);
        add(osmoticBox);
        add(leftMFBox);
        add(rightMFBox);
        
        //panel for all the controls

        getPanel().controlPanel.add(delaySlider.graphic(), vertGBC);
        getPanel().controlPanel.add(temperatureSelect.graphic(), vertGBC);
        getPanel().controlPanel.add(initPanel.graphic(), vertGBC);

        IAction reinitDisplayAction = new IAction() {
        	public void actionPerformed() {
        		molePumpLeft.actionPerformed();
        		leftMFBox.putData(moleFractionAvgLeft.getData());

        		molePumpRight.actionPerformed();
        		rightMFBox.putData(moleFractionAvgRight.getData());
        		osmoticBox.putData(osmosisPMeterAvg.getData());
        		getDisplayBox(sim.box).graphic().repaint();
        	}
        };

        IAction resetDisplayAction = new IAction() {
        	public void actionPerformed() {
        		molePumpLeft.actionPerformed();
        		leftMFBox.putData(moleFractionAvgLeft.getData());

        		molePumpRight.actionPerformed();
        		rightMFBox.putData(moleFractionAvgRight.getData());
        		osmoticBox.putData(osmosisPMeterAvg.getData());
        	}
        };

        getController().getReinitButton().setPostAction(reinitDisplayAction);
        getController().getResetAveragesButton().setPostAction(resetDisplayAction);

        if (sim.getSpace() instanceof Space3D) { // 3D
            ((etomica.graphics.DisplayBoxCanvasG3DSys)displayBox.canvas).setBackgroundColor(Color.WHITE);
            ((etomica.graphics.DisplayBoxCanvasG3DSys)displayBox.canvas).setBoundaryFrameColor(Color.BLACK);
            ((etomica.graphics.DisplayBoxCanvasG3DSys)displayBox.canvas).setPlaneColor(Color.GREEN);
        }

    }

	//drawable that puts a line down the middle of the box, where the
	//semipermeable membrane potential acts
    protected class MyWall implements Drawable {
    	public void draw(Graphics g, int[] origin, double scale) {
    		if(sim.getSpace() instanceof Space2D) {
    		    int x1 = origin[0]+(int)(0.5*scale*sim.box.getBoundary().getBoxSize().getX(0));
    		    int y1 = origin[1];
			    int h = (int)(scale*sim.box.getBoundary().getBoxSize().getX(1));
			    int w = 4;
			    g.setColor(Color.green);
    		    g.fillRect(x1-w, y1, w, h);
    		}
    	}
    }

    public static void main(String[] args) {

        OsmosisSim sim = null;

        Space sp = Space3D.getInstance();
        sim = new OsmosisSim(sp);

        sim.getController().setSleepPeriod(1);
        sim.getController().addActivity(new ActivityIntegrate(sim.integrator));

        Osmosis osmosis = new Osmosis(sim, sp);
        SimulationGraphic.makeAndDisplayFrame(osmosis.getPanel(), APP_NAME);
    }

    protected class InitializeMolecules {

    	private int MINIMUM = 0;
    	private int MOLECULE_MAX = 150;
    	private int SOLUTE_MAX = 50;
    	private int SOLUTE_ON_LEFT_MAX = 100;

    	private JPanel mainPanel;
    	private DeviceSlider total;
    	private DeviceSlider soluteVsSolvent;
    	private DeviceSlider soluteOnLeft;
    	private ConfigurationLatticeWithPlane config;
    	private SimulationRestart simRestart;
    	
    	private int speciesSolventTotal = OsmosisSim.initialSolvent;
    	private int speciesSoluteTotal = OsmosisSim.initialSolute;
    	private float solutePct = (float)speciesSoluteTotal / (float)(speciesSolventTotal + speciesSoluteTotal);

    	public InitializeMolecules(ConfigurationLatticeWithPlane configuration,
    			                   SimulationRestart restart) {

    		this.config = configuration;
    		this.simRestart = restart;

    		// Slider that selects total number of molecules
    		JPanel totalPanel = new JPanel();
            totalPanel.setBorder(new TitledBorder(null, "Total Number of Molecules", TitledBorder.CENTER, TitledBorder.TOP));
    		total = new DeviceSlider(sim.getController());
    		total.setMinimum(MINIMUM);
    		total.setMaximum(MOLECULE_MAX);
    		total.setNMajor(5);
    		total.setShowValues(false);

    		total.setValue(speciesSolventTotal + speciesSoluteTotal);
    		totalPanel.add(total.graphic());

    		ChangeListener totalChange = new ChangeListener() {
    			public void stateChanged(ChangeEvent evt) {

    	    		IAction setAction = new IAction() {
    	    			public void actionPerformed() {
    	    				speciesSoluteTotal = Math.round(((float)total.getValue()) *
    					               (((float)soluteVsSolvent.getValue()) / 100.0f));
    			            speciesSolventTotal = (int)total.getValue() - speciesSoluteTotal;
    			            sim.box.setNMolecules(sim.speciesSolvent, speciesSolventTotal);
    			            sim.box.setNMolecules(sim.speciesSolute, speciesSoluteTotal);
    	    				simRestart.getDataResetAction().actionPerformed();
    	    				simRestart.actionPerformed();
    	            		getDisplayBox(sim.box).graphic().repaint();
    	    			}
    	    		};

    	    		// Need to pause controller, do action, resume controller
    	    		// which is why the action is implemented in this manner.
    	    		sim.getController().submitActionInterrupt(setAction);
    			}
    		};

    		total.getSlider().addChangeListener(totalChange);


    		// Slider that selects solute vx solvent ratio
    		JPanel soluteVsSolventPanel = new JPanel();
    		soluteVsSolventPanel.setBorder(new TitledBorder(null, "Percentage of Solute(vs Solvent)", TitledBorder.CENTER, TitledBorder.TOP));
    		soluteVsSolvent = new DeviceSlider(sim.getController());
    		soluteVsSolvent.setMinimum(MINIMUM);
    		soluteVsSolvent.setMaximum(SOLUTE_MAX);
    		soluteVsSolvent.setNMajor(5);
    		soluteVsSolvent.setShowValues(false);
    		soluteVsSolvent.setValue(Math.round(solutePct * 100.0f));
    		soluteVsSolventPanel.add(soluteVsSolvent.graphic());

            soluteVsSolvent.getSlider().addChangeListener(totalChange);


    		//Slider that selects pct. of solute left of plane
    		JPanel soluteOnLeftPanel = new JPanel();
    		soluteOnLeftPanel.setBorder(new TitledBorder(null, "Percentage of Solute on Left", TitledBorder.CENTER, TitledBorder.TOP));
    		soluteOnLeft = new DeviceSlider(sim.getController());
    		soluteOnLeft.setMinimum(MINIMUM);
    		soluteOnLeft.setMaximum(SOLUTE_ON_LEFT_MAX);
    		soluteOnLeft.setNMajor(4);
    		soluteOnLeft.setShowValues(false);
    		soluteOnLeft.setValue(50);
    		soluteOnLeftPanel.add(soluteOnLeft.graphic());

    		ChangeListener pctChange = new ChangeListener() {
    			public void stateChanged(ChangeEvent evt) {

    	    		IAction setAction = new IAction() {
    	    			public void actionPerformed() {
                            config.setSpeciesAllocation(sim.speciesSolute, (((float)soluteOnLeft.getValue()) / 100.0f));
    				        simRestart.getDataResetAction().actionPerformed();
    				        simRestart.actionPerformed();
    				        getDisplayBox(sim.box).graphic().repaint();
    	    			}
    	    		};

    	    		// Need to pause controller, do action, resume controller
    	    		// which is why the action is implemented in this manner.
    	    		sim.getController().submitActionInterrupt(setAction);
    			}
    		};

    		soluteOnLeft.getSlider().addChangeListener(pctChange);


            mainPanel = new JPanel(new GridLayout(0, 1));
    		mainPanel.add(totalPanel);
    		mainPanel.add(soluteVsSolventPanel);
    		mainPanel.add(soluteOnLeftPanel);
    	}


    	public JPanel graphic() {
    		return mainPanel;
    	}

    	public void addStateChangedListener(ChangeListener cl) {
    		total.getSlider().addChangeListener(cl);
    		soluteVsSolvent.getSlider().addChangeListener(cl);;
        	soluteOnLeft.getSlider().addChangeListener(cl);
    	}

    }
}
