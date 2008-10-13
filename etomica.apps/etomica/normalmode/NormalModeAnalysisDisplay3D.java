package etomica.normalmode;

import etomica.action.activity.ActivityIntegrate;
import etomica.api.IBox;
import etomica.box.Box;
import etomica.graphics.SimulationGraphic;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Pixel;

/**
 * 
 * 
 * @author Tai Boon Tan
 */
public class NormalModeAnalysisDisplay3D extends Simulation {

    private static final long serialVersionUID = 1L;
	private static final String APP_NAME = "3-D Harmonic Oscillator";

	public NormalModeAnalysisDisplay3D(Space _space, int numAtoms, double density, double 
			temperature, String filename){
        super(_space, true);
     
        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        getSpeciesManager().addSpecies(species);
        
        box = new Box(this, space);
        addBox(box);
        box.setNMolecules(species, numAtoms);
        
        double L = Math.pow(4.0/density, 1.0/3.0);
        primitive = new PrimitiveCubic(space, L);
        int n = (int)Math.round(Math.pow(numAtoms/4, 1.0/3.0));
        nCells = new int[]{n, n, n};
        boundary = new BoundaryRectangularPeriodic(space, getRandom(), n*L);
        Basis basis = new BasisCubicFcc();
        box.setBoundary(boundary);

        coordinateDefinition = new CoordinateDefinitionLeaf(this, box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(nCells);
        
        String fileName = "CB_FCC_n12_T01_Mode01";
        normalModes = new NormalModesFromFile(fileName, space.D());
        normalModes.setTemperature(temperature);
                
        WaveVectorFactory waveVectorFactory = normalModes.getWaveVectorFactory();
        waveVectorFactory.makeWaveVectors(box);
                
        integrator = new IntegratorHarmonic(random, 0.00001, temperature, space);

        integrator.setOmegaSquared(normalModes.getOmegaSquared(box), waveVectorFactory.getCoefficients());
        integrator.setEigenVectors(normalModes.getEigenvectors(box));
        integrator.setWaveVectors(waveVectorFactory.getWaveVectors());
        integrator.setWaveVectorCoefficients(waveVectorFactory.getCoefficients());
        integrator.setCoordinateDefinition(coordinateDefinition);
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
 
        int numAtoms =32;
        double density = 1.256;
        double temperature = 0.1;
        //instantiate simulation
        NormalModeAnalysisDisplay3D sim = new NormalModeAnalysisDisplay3D(Space.getInstance(3), 
        		numAtoms, density, temperature, APP_NAME);
        
        SimulationGraphic simGraphic = new SimulationGraphic(sim, Space.getInstance(3), sim.getController());
        simGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(50));
        simGraphic.makeAndDisplayFrame(APP_NAME);
        
	}
	
	public static class Applet extends javax.swing.JApplet{
		
        int numAtoms =108;
        double density = 1.256;
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
	protected NormalModes normalModes;
	protected CoordinateDefinitionLeaf coordinateDefinition;
	
}
