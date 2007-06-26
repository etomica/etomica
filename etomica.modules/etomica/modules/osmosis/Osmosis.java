package etomica.modules.osmosis;

import java.awt.Color;
import java.awt.FlowLayout;
import java.awt.Graphics;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.util.ArrayList;

import javax.swing.JPanel;
import javax.swing.border.TitledBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import etomica.action.Action;
import etomica.action.SimulationRestart;
import etomica.config.ConfigurationLatticeWithPlane;
import etomica.data.AccumulatorAverage;
import etomica.data.DataPump;
import etomica.data.DataSourceCountTime;
import etomica.data.meter.MeterLocalMoleFraction;
import etomica.data.meter.MeterTemperature;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DeviceThermoSelector;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxesCAE;
import etomica.graphics.DisplayPhase;
import etomica.graphics.DisplayTimer;
import etomica.graphics.Drawable;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.graphics.DisplayBox.LabelType;
import etomica.lattice.LatticeCubicSimple;
import etomica.math.geometry.Cuboid;
import etomica.math.geometry.Plane;
import etomica.math.geometry.Rectangle;
import etomica.potential.P1HardBoundary;
import etomica.space.IVector;
import etomica.space2d.Space2D;
import etomica.space2d.Vector2D;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.units.Kelvin;
import etomica.units.Unit;
import etomica.util.Constants.CompassDirection;

/**
 * Osmosis module.
 * @author Jhumpa Adhikari
 * @author Andrew Schultz
 */

public class Osmosis extends SimulationGraphic {

	private final static String APP_NAME = "Osmosis";
	private final static int REPAINT_INTERVAL = 40;

    public DataSourceCountTime cycles;
    public DisplayBox displayCycles;
    public MeterOsmoticPressure osmosisPMeter;
    public MeterLocalMoleFraction moleFractionRight, moleFractionLeft;
    public OsmosisSim sim;
    
    public Osmosis(OsmosisSim simulation) {

    	super(simulation, GRAPHIC_ONLY, APP_NAME, REPAINT_INTERVAL);

        ArrayList dataStreamPumps = getController().getDataStreamPumps();
        
    	sim = simulation;
        final int thickness = 4;

        Unit tUnit = Kelvin.UNIT;

        ConfigurationLatticeWithPlane config = null;

	    //display of phase
        final DisplayPhase displayPhase = getDisplayPhase(sim.phase);
        ColorSchemeByType colorScheme = new ColorSchemeByType();

        colorScheme.setColor(sim.speciesSolvent.getMoleculeType(), Color.blue);
        colorScheme.setColor(sim.speciesSolute.getMoleculeType(), Color.red);
        displayPhase.setColorScheme(colorScheme);
        displayPhase.setAlign(1,DisplayPhase.CENTER);
        displayPhase.setOriginShift(0, thickness);
        displayPhase.setOriginShift(1, -thickness);
        if (sim.getSpace() instanceof Space2D) {
            displayPhase.addDrawable(new MyWall());
            config = new ConfigurationLatticeWithPlane(new LatticeCubicSimple(2, 1.0), null);
        }
        else if (sim.getSpace() instanceof Space3D) {
        	Plane plane = new Plane(sim.getSpace());
        	((etomica.graphics.DisplayPhaseCanvasG3DSys)displayPhase.canvas).addPlane(plane);
            config = new ConfigurationLatticeWithPlane(new LatticeCubicSimple(3, 1.0), plane); 
            config.addSpecies((etomica.species.Species)sim.speciesSolvent);
            config.addSpecies((etomica.species.Species)sim.speciesSolute);
        }

        config.initializeCoordinates(sim.phase);

        SimulationRestart simRestart = getController().getSimRestart();
        
        simRestart.setConfiguration(config);
        getController().setShape("VERTICAL"); //three choices "HORIZONTAL","AUTOMATIC"           

        cycles = new DataSourceCountTime();
        displayCycles = new DisplayTimer(sim.integrator);
        displayCycles.setLabelType(LabelType.BORDER);
        displayCycles.setLabel("Cycle Time");
	    displayCycles.setPrecision(6);	

	    // Right side of membrane osmotic
        osmosisPMeter = new MeterOsmoticPressure(sim.getSpace(), new P1HardBoundary[]{sim.boundaryHardA,sim.boundaryHardB}); 
        osmosisPMeter.setIntegrator(sim.integrator);
        final AccumulatorAverage osmosisPMeterAvg = new AccumulatorAverage();
        final DataPump osmosisPump = new DataPump(osmosisPMeter, osmosisPMeterAvg);
        dataStreamPumps.add(osmosisPump);
        sim.integrator.addIntervalAction(osmosisPump);
        sim.integrator.setActionInterval(osmosisPump, 40);
        final DisplayBoxesCAE osmoticBox = new DisplayBoxesCAE();
        osmoticBox.setAccumulator(osmosisPMeterAvg);
        osmoticBox.setPrecision(6);

        //
        // temperature panel
        //

	    DeviceThermoSelector tSelect = new DeviceThermoSelector(sim, sim.integrator);
	    tSelect.setTemperatures(new double[] {50.,100.,300.,600.,1000.});
	    tSelect.setUnit(tUnit);
	    tSelect.setSelected(0); //sets adiabatic as selected temperature
		MeterTemperature thermometer = new MeterTemperature();
		thermometer.setPhase(sim.phase);
		DisplayBox tBox = new DisplayBox();
        DataPump tempPump = new DataPump(thermometer, tBox);
        sim.integrator.addIntervalAction(tempPump);
		tBox.setUnit(tUnit);
		tBox.setLabel("Measured value");
		tBox.setLabelPosition(CompassDirection.NORTH);
	    tSelect.getLabel().setText("Set value");

        JPanel temperaturePanel = new JPanel(new GridBagLayout());

        temperaturePanel.setBorder(new TitledBorder(null, "Temperature (K)", TitledBorder.CENTER, TitledBorder.TOP));
        temperaturePanel.add(tSelect.graphic(null),SimulationPanel.getHorizGBC());
        temperaturePanel.add(tBox.graphic(null),SimulationPanel.getHorizGBC());

        // Right side of membrane mole fraction
        moleFractionRight = new MeterLocalMoleFraction();
        moleFractionRight.setPhase(sim.phase);
        IVector dimensions = sim.phase.getBoundary().getDimensions();

        if (sim.getSpace() instanceof Space2D) { // 2D
            moleFractionRight.setShape(new Rectangle(sim.getSpace(), dimensions.x(0)*0.5, dimensions.x(1)));
            moleFractionRight.setShapeOrigin(new Vector2D(dimensions.x(0)*0.25, 0));
        }
        else if (sim.getSpace() instanceof Space3D) { // 3D
            moleFractionRight.setShape(new Cuboid(sim.getSpace(), dimensions.x(0)*0.5, dimensions.x(1), dimensions.x(2)));
            moleFractionRight.setShapeOrigin(new Vector3D(dimensions.x(0)*0.25, 0, 0));
        }

        moleFractionRight.setSpecies(sim.speciesSolute);
        final AccumulatorAverage moleFractionAvgRight = new AccumulatorAverage();
        final DataPump molePumpRight = new DataPump(moleFractionRight, moleFractionAvgRight);
        dataStreamPumps.add(molePumpRight);
        sim.integrator.addIntervalAction(molePumpRight);
        final DisplayBoxesCAE rightMFBox = new DisplayBoxesCAE();
        rightMFBox.setAccumulator(moleFractionAvgRight);
        rightMFBox.setPrecision(8);

        // Left side of membrane mole fraction
        moleFractionLeft = new MeterLocalMoleFraction();
        moleFractionLeft.setPhase(sim.phase);

        if (sim.getSpace() instanceof Space2D) { // 2D
            moleFractionLeft.setShape(new Rectangle(sim.getSpace(), dimensions.x(0)*0.5, dimensions.x(1)));
            moleFractionLeft.setShapeOrigin(new Vector2D(-dimensions.x(0)*0.25, 0));
        }
        else if (sim.getSpace() instanceof Space3D) { // 3D
            moleFractionLeft.setShape(new Cuboid(sim.getSpace(), dimensions.x(0)*0.5, dimensions.x(1), dimensions.x(2)));

            moleFractionLeft.setShapeOrigin(new Vector3D(-dimensions.x(0)*0.25, 0, 0));
        }

        moleFractionLeft.setSpecies(sim.speciesSolute);
        final AccumulatorAverage moleFractionAvgLeft = new AccumulatorAverage();
        final DataPump molePumpLeft = new DataPump(moleFractionLeft, moleFractionAvgLeft);
        dataStreamPumps.add(molePumpLeft);
        sim.integrator.addIntervalAction(molePumpLeft);
        final DisplayBoxesCAE leftMFBox = new DisplayBoxesCAE();
        leftMFBox.setAccumulator(moleFractionAvgLeft);
        leftMFBox.setPrecision(8);

        ChangeListener cl = new ChangeListener() {
        	public void stateChanged(ChangeEvent evt) {
				molePumpLeft.actionPerformed();
				leftMFBox.putData(moleFractionAvgLeft.getData());
				leftMFBox.repaint();
		
				molePumpRight.actionPerformed();
				rightMFBox.putData(moleFractionAvgRight.getData());
				rightMFBox.repaint();
				osmoticBox.putData(osmosisPMeterAvg.getData());
				osmoticBox.repaint();
				getDisplayPhase(sim.phase).graphic().repaint();
        	}
        };

        InitializeMolecules initPanel = new InitializeMolecules(config, simRestart);
        initPanel.addStateChangedListener(cl);

        // panel for osmotic pressure
        JPanel osmoticPanel = new JPanel(new FlowLayout());
        osmoticPanel.setBorder(new TitledBorder(null, "Osmotic Pressure (PV/Nk)", TitledBorder.CENTER, TitledBorder.TOP));
        osmoticPanel.add(osmoticBox.graphic(null));

        // left side panel for mole fraction
        JPanel leftMoleFractionPanel = new JPanel(new FlowLayout());
        leftMoleFractionPanel.setBorder(new TitledBorder(null, "Mole Fraction (nSolute/nSolution)", TitledBorder.CENTER, TitledBorder.TOP));
        leftMoleFractionPanel.add(leftMFBox.graphic(null));

        // left side metrics panel
        JPanel leftMetricsPanel = new JPanel(new GridLayout(0, 1));
        leftMetricsPanel.setBorder(new TitledBorder(null, "Left of Membrane", TitledBorder.CENTER, TitledBorder.TOP));
        leftMetricsPanel.add(leftMoleFractionPanel);

        // right side panel for mole fraction
        JPanel rightMoleFractionPanel = new JPanel(new FlowLayout());
        rightMoleFractionPanel.setBorder(new TitledBorder(null, "Mole Fraction (nSolute/nSolution)", TitledBorder.CENTER, TitledBorder.TOP));
        rightMoleFractionPanel.add(rightMFBox.graphic(null));

        // right side metrics panel
        JPanel rightMetricsPanel = new JPanel(new GridLayout(0, 1));
        rightMetricsPanel.setBorder(new TitledBorder(null, "Right of Membrane", TitledBorder.CENTER, TitledBorder.TOP));
        rightMetricsPanel.add(rightMoleFractionPanel);


        //panel for all the controls

        GridBagConstraints vertGBC = SimulationPanel.getVertGBC();

        getPanel().controlPanel.add(temperaturePanel, vertGBC);
        getPanel().controlPanel.add(initPanel.graphic(), vertGBC);
        getPanel().plotPanel.add(displayCycles.graphic(), vertGBC);
        getPanel().plotPanel.add(osmoticPanel, vertGBC);
        getPanel().plotPanel.add(leftMetricsPanel, vertGBC);
        getPanel().plotPanel.add(rightMetricsPanel, vertGBC);

        Action reinitDisplayAction = new Action() {
        	public void actionPerformed() {
        		molePumpLeft.actionPerformed();
        		leftMFBox.putData(moleFractionAvgLeft.getData());
        		leftMFBox.repaint();

        		molePumpRight.actionPerformed();
        		rightMFBox.putData(moleFractionAvgRight.getData());
        		rightMFBox.repaint();
        		osmoticBox.putData(osmosisPMeterAvg.getData());
        		osmoticBox.repaint();
        		getDisplayPhase(sim.phase).graphic().repaint();
        	}
        };

        Action resetDisplayAction = new Action() {
        	public void actionPerformed() {
        		molePumpLeft.actionPerformed();
        		leftMFBox.putData(moleFractionAvgLeft.getData());
        		leftMFBox.repaint();

        		molePumpRight.actionPerformed();
        		rightMFBox.putData(moleFractionAvgRight.getData());
        		rightMFBox.repaint();
        		osmoticBox.putData(osmosisPMeterAvg.getData());
        		osmoticBox.repaint();
        	}
        };

        getController().getReinitButton().setPostAction(reinitDisplayAction);
        getController().getResetAveragesButton().setPostAction(resetDisplayAction);

        if (sim.getSpace() instanceof Space3D) { // 3D
            ((etomica.graphics.DisplayPhaseCanvasG3DSys)displayPhase.canvas).setBackgroundColor(Color.WHITE);
            ((etomica.graphics.DisplayPhaseCanvasG3DSys)displayPhase.canvas).setBoundaryFrameColor(Color.BLACK);
            ((etomica.graphics.DisplayPhaseCanvasG3DSys)displayPhase.canvas).setPlaneColor(Color.GREEN);
        }

    }

	//drawable that puts a line down the middle of the box, where the
	//semipermeable membrane potential acts
    protected class MyWall implements Drawable {
    	public void draw(Graphics g, int[] origin, double scale) {
    		if(sim.getSpace() instanceof Space2D) {
    		    int x1 = origin[0]+(int)(0.5*scale*sim.phase.getBoundary().getDimensions().x(0));
    		    int y1 = origin[1];
			    int h = (int)(scale*sim.phase.getBoundary().getDimensions().x(1));
			    int w = 4;
			    g.setColor(Color.green);
    		    g.fillRect(x1-w, y1, w, h);
    		}
    	}
    }

    public static void main(String[] args) {

        OsmosisSim sim = null;

        if(true) { // 3D Case
    	    sim = new OsmosisSim(Space3D.getInstance());
        }
        else { // 2D case
        	sim = new OsmosisSim(Space2D.getInstance());
        }

        sim.activityIntegrate.setDoSleep(true);
        sim.activityIntegrate.setSleepPeriod(1);

        Osmosis osmosis = new Osmosis(sim);
        SimulationGraphic.makeAndDisplayFrame(osmosis.getPanel(), APP_NAME);
    }

    public static class Applet extends javax.swing.JApplet {
	    public void init() {

	    	OsmosisSim sim = null;

	        if(false) { // 3D Case
	    	    sim = new OsmosisSim(Space3D.getInstance());
	        }
	        else { // 2D case
	        	sim = new OsmosisSim(Space2D.getInstance());
	        }

	        sim.activityIntegrate.setDoSleep(true);
	        sim.activityIntegrate.setSleepPeriod(1);

		    getContentPane().add(new Osmosis(sim).getPanel());
	    }

        private static final long serialVersionUID = 1L;
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
    				speciesSoluteTotal = Math.round(((float)total.getValue()) *
    						               (((float)soluteVsSolvent.getValue()) / 100.0f));
    				speciesSolventTotal = (int)total.getValue() - speciesSoluteTotal;
    				
    				sim.phase.getAgent(sim.speciesSolvent).setNMolecules(speciesSolventTotal);
    				sim.phase.getAgent(sim.speciesSolute).setNMolecules(speciesSoluteTotal);
    				simRestart.getDataResetAction().actionPerformed();
    				simRestart.actionPerformed();
    				getDisplayPhase(sim.phase).graphic().repaint();
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
    		soluteOnLeftPanel.add(soluteOnLeft.graphic());

    		ChangeListener pctChange = new ChangeListener() {
    			public void stateChanged(ChangeEvent evt) {
    				config.setSpeciesAllocation(sim.speciesSolute, (((float)soluteOnLeft.getValue()) / 100.0f));
    				simRestart.getDataResetAction().actionPerformed();
    				simRestart.actionPerformed();
    				getDisplayPhase(sim.phase).graphic().repaint();
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
