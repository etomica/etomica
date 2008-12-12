package etomica.normalmode;

import etomica.action.activity.ActivityIntegrate;
import etomica.api.IBox;
import etomica.box.Box;
import etomica.graphics.SimulationGraphic;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;

/**
 * 
 * 
 * @author Tai Boon Tan
 */
public class NormalModeAnalysisDisplay1D extends Simulation {

    private static final long serialVersionUID = 1L;
	private static final String APP_NAME = "1-D Harmonic Oscillator";

	public NormalModeAnalysisDisplay1D(Space space){
        super(space, true);
     
        species = new SpeciesSpheresMono(this, space);
        getSpeciesManager().addSpecies(species);
        
        box = new Box(this, space);
        addBox(box);
        box.setNMolecules(species, numAtoms);
        
        primitive = new PrimitiveCubic(space, 1.0/density);
        boundary = new BoundaryRectangularPeriodic(space, getRandom(), numAtoms/density);
        nCells = new int[]{numAtoms};
        box.setBoundary(boundary);

        coordinateDefinition = new CoordinateDefinitionLeaf(this, box, primitive, space);
        coordinateDefinition.initializeCoordinates(nCells);
        
        nm = new NormalModes1DHR(space.D());
        nm.setTemperature(temperature);
        
        waveVectorFactory = nm.getWaveVectorFactory();
        waveVectorFactory.makeWaveVectors(box);
        
        integrator = new IntegratorHarmonic(random, 0.01, temperature, space);
        integrator.setCoordinateDefinition(coordinateDefinition);
        integrator.setWaveVectors(waveVectorFactory.getWaveVectors());
        integrator.setWaveVectorCoefficients(waveVectorFactory.getCoefficients());
        integrator.setOmegaSquared(nm.getOmegaSquared(box), waveVectorFactory.getCoefficients());
        integrator.setEigenVectors(nm.getEigenvectors(box));
        integrator.setTemperature(temperature);
        
        
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        activityIntegrate.setSleepPeriod(1);
        
        getController().addAction(activityIntegrate);
        integrator.setBox(box);
        
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
 
        //instantiate simulation
		Space sp = Space.getInstance(1);
        NormalModeAnalysisDisplay1D sim = new NormalModeAnalysisDisplay1D(sp);
        
        SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, sp, sim.getController());
        simGraphic.makeAndDisplayFrame(APP_NAME);
        
	}
	
	public static class Applet extends javax.swing.JApplet{
		
        int numAtoms =30;
        double density = 1.0;
        double temperature = 0.1;
		
		public void init(){
	
		}
	}
	
	
	
	protected IntegratorHarmonic integrator;
	protected ActivityIntegrate activityIntegrate;
	protected IBox box;
	protected Boundary boundary;
	protected Primitive primitive;
	protected Basis basis;
	protected int[] nCells;
	protected SpeciesSpheresMono species;
	protected NormalModes1DHR nm;
	protected CoordinateDefinitionLeaf coordinateDefinition;
	protected WaveVectorFactory waveVectorFactory;
	//Default:
	protected double density = 0.5;
	protected double temperature = 1.0;
	protected int numAtoms = 50;




	
}
