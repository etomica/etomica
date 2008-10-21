package etomica.normalmode;

import etomica.action.activity.ActivityIntegrate;
import etomica.api.IBox;
import etomica.api.IVector;
import etomica.box.Box;
import etomica.graphics.SimulationGraphic;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.BasisOrthorhombicHexagonal;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveOrthorhombicHexagonal;
import etomica.paracetamol.BasisOrthorhombicParacetamol;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space2d.Space2D;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Pixel;

/**
 * 
 * Harmonic anaylsis for 2-D soft-sphere model
 * 
 * @author Tai Boon Tan
 */
public class NormalModeAnalysisDisplay2D extends Simulation {

    private static final long serialVersionUID = 1L;
	private static final String APP_NAME = "2-D Harmonic Oscillator";

	public NormalModeAnalysisDisplay2D(Space _space, int numAtoms, int[] nCells, double 
			temperature, String filename){
        super(_space, true);
     
        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        getSpeciesManager().addSpecies(species);
        
        box = new Box(this, space);
        addBox(box);
        box.setNMolecules(species, numAtoms);
        
        primitive = new PrimitiveOrthorhombicHexagonal(space, 1);
        IVector[] dimension = space.makeVectorArray(space.D());
        for (int i=0; i<space.D(); i++){
        	dimension[i].Ea1Tv1(nCells[i], primitive.vectors()[i]);
        }
        boundary = new BoundaryDeformablePeriodic(space, random, dimension);
        Basis basis = new BasisOrthorhombicHexagonal();
        box.setBoundary(boundary);

        coordinateDefinition = new CoordinateDefinitionLeaf(this, box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(nCells);
        
        String fileName = "2D_CB_FCC_n12_T01";
        normalModes = new NormalModesFromFile(fileName, space.D());
        normalModes.setTemperature(temperature);
                
        WaveVectorFactory waveVectorFactory = normalModes.getWaveVectorFactory();
        waveVectorFactory.makeWaveVectors(box);
                
        integrator = new IntegratorHarmonic(random, 0.000001, temperature, space);

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
 
        double temperature = 0.1;
        Space sp = Space2D.getInstance();
        int[] dimension = new int[] {2,1};
        int numAtoms = 2*dimension[0]*dimension[1];
        //instantiate simulation
        NormalModeAnalysisDisplay2D sim = new NormalModeAnalysisDisplay2D(sp, 
        		numAtoms, dimension, temperature, APP_NAME);
        
        SimulationGraphic simGraphic = new SimulationGraphic(sim, sim.space, sim.getController());
        simGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(30));
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
