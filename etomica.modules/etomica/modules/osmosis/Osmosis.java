package etomica.modules.osmosis;

import java.awt.Color;
import java.awt.FlowLayout;
import java.awt.Graphics;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;

import javax.swing.JPanel;
import javax.swing.border.TitledBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import etomica.action.Action;
import etomica.action.SimulationRestart;
import etomica.atom.AtomTypeSphere;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverage;
import etomica.data.DataPump;
import etomica.data.DataSourceCountTime;
import etomica.data.meter.MeterLocalMoleFraction;
import etomica.data.meter.MeterTemperature;
import etomica.exception.ConfigurationOverlapException;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DeviceNSelector;
import etomica.graphics.DeviceThermoSelector;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxesCAE;
import etomica.graphics.DisplayPhase;
import etomica.graphics.DisplayTimer;
import etomica.graphics.Drawable;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.graphics.DisplayBox.LabelType;
import etomica.integrator.Integrator;
import etomica.lattice.LatticeCubicSimple;
import etomica.math.geometry.Cuboid;
import etomica.math.geometry.Rectangle;
import etomica.modifier.Modifier;
import etomica.potential.P1HardBoundary;
import etomica.potential.P2HardSphere;
import etomica.space.IVector;
import etomica.space2d.Space2D;
import etomica.space3d.Space3D;
import etomica.space2d.Vector2D;
import etomica.space3d.Vector3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Dimension;
import etomica.units.Kelvin;
import etomica.units.Length;
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
    public MeterLocalMoleFraction moleFraction;
    public OsmosisSim sim;

    public Osmosis(OsmosisSim simulation) {

    	super(simulation, GRAPHIC_ONLY, APP_NAME, REPAINT_INTERVAL);

    	sim = simulation;

        final int thickness = 4;

        Unit tUnit = Kelvin.UNIT;

        ConfigurationLattice config = new ConfigurationLattice(new LatticeCubicSimple(2, 1.0));

        SimulationRestart simRestart = getController().getSimRestart();
        
        simRestart.setConfiguration(config);
        getController().setShape("VERTICAL"); //three choices "HORIZONTAL","AUTOMATIC"           

	    //display of phase
        final DisplayPhase displayPhase = getDisplayPhase(sim.phase);
        ColorSchemeByType colorScheme = new ColorSchemeByType();

        colorScheme.setColor(sim.speciesA.getMoleculeType(), Color.blue);
        colorScheme.setColor(sim.speciesB.getMoleculeType(), Color.red);
        displayPhase.setColorScheme(colorScheme);
        displayPhase.setAlign(1,DisplayPhase.CENTER);
        displayPhase.setOriginShift(0, thickness);
        displayPhase.setOriginShift(1, -thickness);
        if (sim.getSpace() instanceof Space2D) {
            displayPhase.addDrawable(new MyWall());
        }
        else if (sim.getSpace() instanceof Space3D) {
        	etomica.math.geometry.Plane plane = new etomica.math.geometry.Plane(sim.getSpace());
        	((etomica.graphics.DisplayPhaseCanvasG3DSys)displayPhase.canvas).addPlane(plane);
        }

        cycles = new DataSourceCountTime();
        displayCycles = new DisplayTimer(sim.integrator);
        displayCycles.setLabelType(LabelType.BORDER);
        displayCycles.setLabel("Cycle Time");
	    displayCycles.setPrecision(6);	

        osmosisPMeter = new MeterOsmoticPressure(sim.getSpace(), new P1HardBoundary[]{sim.boundaryHardLeftA}, 
                new P1HardBoundary[]{sim.boundaryHardRightA, sim.boundaryHardB});
        osmosisPMeter.setIntegrator(sim.integrator);
        final AccumulatorAverage osmosisPMeterAvg = new AccumulatorAverage();
        DataPump pump = new DataPump(osmosisPMeter, osmosisPMeterAvg);
        sim.register(osmosisPMeter, pump);
        sim.integrator.addIntervalAction(pump);
        sim.integrator.setActionInterval(pump, 40);
        final DisplayBoxesCAE dBox = new DisplayBoxesCAE();
        dBox.setAccumulator(osmosisPMeterAvg);
        dBox.setPrecision(6);

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
        pump = new DataPump(thermometer, tBox);
        sim.integrator.addIntervalAction(pump);
		tBox.setUnit(tUnit);
		tBox.setLabel("Measured value");
		tBox.setLabelPosition(CompassDirection.NORTH);
	    tSelect.getLabel().setText("Set value");

        JPanel temperaturePanel = new JPanel(new GridBagLayout());

        temperaturePanel.setBorder(new TitledBorder(null, "Temperature (K)", TitledBorder.CENTER, TitledBorder.TOP));
        temperaturePanel.add(tSelect.graphic(null),SimulationPanel.getHorizGBC());
        temperaturePanel.add(tBox.graphic(null),SimulationPanel.getHorizGBC());

        moleFraction = new MeterLocalMoleFraction();
        moleFraction.setPhase(sim.phase);
        IVector dimensions = sim.phase.getBoundary().getDimensions();

        if (sim.getSpace() instanceof Space2D) { // 2D
            moleFraction.setShape(new Rectangle(sim.getSpace(), dimensions.x(0)*0.5, dimensions.x(1)));
            moleFraction.setShapeOrigin(new Vector2D(dimensions.x(0)*0.25, 0));
        }
        else if (sim.getSpace() instanceof Space3D) { // 3D
            moleFraction.setShape(new Cuboid(sim.getSpace(), dimensions.x(0)*0.5, dimensions.x(1), dimensions.x(2)));
            moleFraction.setShapeOrigin(new Vector3D(dimensions.x(0)*0.25, 0, 0));
        }

        moleFraction.setSpecies(sim.speciesB);
        final AccumulatorAverage moleFractionAvg = new AccumulatorAverage();
        pump = new DataPump(moleFraction, moleFractionAvg);
        sim.register(moleFraction, pump);
        sim.integrator.addIntervalAction(pump);
        final DisplayBoxesCAE mfBox = new DisplayBoxesCAE();
        mfBox.setAccumulator(moleFractionAvg);
        mfBox.setPrecision(8);

        DeviceNSelector nASelector = new DeviceNSelector(sim.getController());
        nASelector.setResetAction(simRestart);
        nASelector.setSpeciesAgent(sim.phase.getAgent(sim.speciesA));
        nASelector.getSlider().addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent evt) {
            	displayPhase.repaint();
            }
        });
        nASelector.setMaximum(30);
        
        DeviceNSelector nBSelector = new DeviceNSelector(sim.getController());
        nBSelector.setResetAction(simRestart);
        nBSelector.setSpeciesAgent(sim.phase.getAgent(sim.speciesB));
        nBSelector.getSlider().addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent evt) {
            	displayPhase.repaint();
            }
        });
        nBSelector.setMaximum(10);
	    
        DiameterModifier diameterModifier = new DiameterModifier(sim.potentialAA,sim.potentialBB,sim.potentialAB,
                                                            sim.boundarySemiB,
                                                            sim.boundaryHardTopBottomA,sim.boundaryHardLeftA,sim.boundaryHardRightA,
                                                            sim.boundaryHardB,
                                                            sim.speciesA,sim.speciesB);
        // panel for osmotic pressure
        JPanel osmoticPanel = new JPanel(new FlowLayout());
        osmoticPanel.setBorder(new TitledBorder(null, "Osmotic Pressure (PV/Nk)", TitledBorder.CENTER, TitledBorder.TOP));
        osmoticPanel.add(dBox.graphic(null));

        // panel for mole fraction
        JPanel moleFractionPanel = new JPanel(new FlowLayout());
        TitledBorder titleBorder = new TitledBorder(null, "Mole Fraction (nSolute/nSolution)", TitledBorder.CENTER, TitledBorder.TOP);
        moleFractionPanel.setBorder(titleBorder);
        moleFractionPanel.add(mfBox.graphic(null));

        diameterModifier.setValue(1.0);

		// Solvent molecules slider
        JPanel sliderPanelA = new JPanel(new GridLayout(0,1));
        nASelector.setShowBorder(false);
        sliderPanelA.add(nASelector.graphic(null));
        sliderPanelA.setBorder(new TitledBorder
           (null, "Set "+nASelector.getLabel(), TitledBorder.CENTER, TitledBorder.TOP));

        // Solute molecules slider
        JPanel sliderPanelB = new JPanel(new GridLayout(0,1));
        nBSelector.setShowBorder(false);
        sliderPanelB.add(nBSelector.graphic(null));
        sliderPanelB.setBorder(new TitledBorder
           (null, "Set "+nBSelector.getLabel(), TitledBorder.CENTER, TitledBorder.TOP));

        //panel for all the controls

        GridBagConstraints vertGBC = SimulationPanel.getVertGBC();

        getPanel().controlPanel.add(temperaturePanel, vertGBC);
        getPanel().controlPanel.add(sliderPanelA, vertGBC);
        getPanel().controlPanel.add(sliderPanelB, vertGBC);
        getPanel().plotPanel.add(displayCycles.graphic(), vertGBC);
        getPanel().plotPanel.add(osmoticPanel, vertGBC);
        getPanel().plotPanel.add(moleFractionPanel, vertGBC);

        Action reinitDisplayAction = new Action() {
        	public void actionPerformed() {
        		mfBox.putData(moleFractionAvg.getData());
        		mfBox.repaint();
        		dBox.putData(osmosisPMeterAvg.getData());
        		dBox.repaint();
        		getDisplayPhase(sim.phase).graphic().repaint();
        	}
        };

        Action resetDisplayAction = new Action() {
        	public void actionPerformed() {
        		mfBox.putData(moleFractionAvg.getData());
        		mfBox.repaint();
        		dBox.putData(osmosisPMeterAvg.getData());
        		dBox.repaint();
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

        if(false) { // 3D Case
    	    sim = new OsmosisSim(Space3D.getInstance());
        }
        else { // 2D case
        	sim = new OsmosisSim(Space2D.getInstance());
        }

        sim.activityIntegrate.setDoSleep(true);
        sim.activityIntegrate.setSleepPeriod(1);
        sim.register(sim.integrator);

        Osmosis osmosis = new Osmosis(sim);
        SimulationGraphic.makeAndDisplayFrame(osmosis.getPanel(), APP_NAME);
    }

    
    public class DiameterModifier implements Modifier {
        
        P2HardSphere potentialAA,potentialBB,potentialAB;
        P1HardWall membraneB;
        P1HardBoundary boundaryHard1A,boundaryHard2A,boundaryHard3A;
        P1HardBoundary boundaryHardB;
        SpeciesSpheresMono speciesA,speciesB;
        DisplayPhase display;
        Integrator integrator;
        
        public DiameterModifier(P2HardSphere potentialAA ,P2HardSphere potentialBB ,P2HardSphere potentialAB ,
                          P1HardWall membraneB,
                          P1HardBoundary boundaryHard1A, P1HardBoundary boundaryHard2A, P1HardBoundary boundaryHard3A, 
                          P1HardBoundary boundaryHardB,
                          SpeciesSpheresMono speciesA ,SpeciesSpheresMono speciesB) {
            this.potentialAA = potentialAA;
            this.potentialBB = potentialBB;
            this.potentialAB = potentialAB;
            this.membraneB = membraneB;
            this.boundaryHard1A = boundaryHard1A;
            this.boundaryHard2A = boundaryHard2A;
            this.boundaryHard3A = boundaryHard3A;
            this.boundaryHardB = boundaryHardB;
            this.speciesA = speciesA;
            this.speciesB = speciesB;
        }
        
        public String getLabel() {
            return "a label";
        }
        
        public Dimension getDimension() {
            return Length.DIMENSION;
        }
        
        public void setValue(double d) {
            potentialAA.setCollisionDiameter(d);
            potentialBB.setCollisionDiameter(d);
            potentialAB.setCollisionDiameter(d);
            boundaryHard1A.setCollisionRadius(0.5*d);
            boundaryHard2A.setCollisionRadius(0.5*d);
            boundaryHard3A.setCollisionRadius(0.5*d);
            membraneB.setCollisionRadius(0.5*d);
            boundaryHardB.setCollisionRadius(0.5*d);
            ((AtomTypeSphere)speciesA.getMoleculeType()).setDiameter(d);
            ((AtomTypeSphere)speciesB.getMoleculeType()).setDiameter(d);
            if (display != null) {
                display.repaint();
            }
            if (integrator != null) {
                try {
                    integrator.reset();
                }
                catch (ConfigurationOverlapException e) {
                    //overlaps are likely after increasing the diameter
                }
            }
        }
        
        public double getValue() {
            return ((AtomTypeSphere)speciesA.getMoleculeType()).getDiameter();
        }
        
        public void setDisplay(DisplayPhase newDisplay){
            display = newDisplay;
        }
        
        public void setIntegrator(Integrator newIntegrator) {
            integrator = newIntegrator;
        }
    }

    public static class Applet extends javax.swing.JApplet {
	    public void init() {

	    	OsmosisSim sim = null;

	        if(true) { // 3D Case
	    	    sim = new OsmosisSim(Space3D.getInstance());
	        }
	        else { // 2D case
	        	sim = new OsmosisSim(Space2D.getInstance());
	        }

	        sim.activityIntegrate.setDoSleep(true);
	        sim.activityIntegrate.setSleepPeriod(1);
	        sim.register(sim.integrator);

		    getContentPane().add(new Osmosis(sim).getPanel());
	    }

        private static final long serialVersionUID = 1L;
    }

}








