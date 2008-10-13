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

	public NormalModeAnalysisDisplay1D(Space _space, int numAtoms, double density, double 
			temperature, String filename){
        super(_space, true);
     
        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
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
        
        WaveVectorFactory waveVectorFactory = nm.getWaveVectorFactory();
        waveVectorFactory.makeWaveVectors(box);
        
        integrator = new IntegratorHarmonic(random, 0.001, temperature, space);
        integrator.setCoordinateDefinition(coordinateDefinition);
        integrator.setWaveVectors(waveVectorFactory.getWaveVectors());
        integrator.setWaveVectorCoefficients(waveVectorFactory.getCoefficients());
        integrator.setOmegaSquared(nm.getOmegaSquared(box), waveVectorFactory.getCoefficients());
        integrator.setEigenVectors(nm.getEigenvectors(box));
        integrator.setTemperature(temperature);
        
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        activityIntegrate.setSleepPeriod(0);
        
        getController().addAction(activityIntegrate);
        integrator.setBox(box);
        
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
 
        int numAtoms =20;
        double density = 0.4;
        double temperature = 0.1;
        //instantiate simulation
        NormalModeAnalysisDisplay1D sim = new NormalModeAnalysisDisplay1D(Space.getInstance(1), 
        		numAtoms, density, temperature, APP_NAME);
        
        SimulationGraphic simGraphic = new SimulationGraphic(sim, Space.getInstance(1), sim.getController());
        simGraphic.makeAndDisplayFrame(APP_NAME);
        
	}
	
	public static class Applet extends javax.swing.JApplet{
		
        int numAtoms =20;
        double density = 0.4;
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
	
}
