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
import etomica.space.ISpace;
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

	public NormalModeAnalysisDisplay3D(Space _space){
        super(_space, true);
        this.space = _space;
        
        species = new SpeciesSpheresMono(this, space);
        getSpeciesManager().addSpecies(species);
        
        box = new Box(space);
        addBox(box);
        box.setNMolecules(species, numAtoms);
        
        L = Math.pow(4.0/density, 1.0/3.0);
        primitive = new PrimitiveCubic(space, L);
        n = (int)Math.round(Math.pow(numAtoms/4, 1.0/3.0));
        nCells = new int[]{n, n, n};
        boundary = new BoundaryRectangularPeriodic(space, n*L);
        basis = new BasisCubicFcc();
        box.setBoundary(boundary);

        coordinateDefinition = new CoordinateDefinitionLeaf(this, box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(nCells);
        
        //String fileName = "CB_FCC_n12_T01_Mode01";
        nm = new NormalModes3D(space, primitive, basis);
        nm.setTemperature(temperature);
                
        waveVectorFactory = nm.getWaveVectorFactory();
        waveVectorFactory.makeWaveVectors(box);
                
        integrator = new IntegratorHarmonic(random, 0.0001, temperature, space);

        integrator.setOmegaSquared(nm.getOmegaSquared(box), waveVectorFactory.getCoefficients());
        integrator.setEigenVectors(nm.getEigenvectors(box));
        integrator.setWaveVectors(waveVectorFactory.getWaveVectors());
        integrator.setWaveVectorCoefficients(waveVectorFactory.getCoefficients());
        integrator.setCoordinateDefinition(coordinateDefinition);
        integrator.setTemperature(temperature);
        integrator.setOneWV(true);
        integrator.setWaveVectorNum(0);
        
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        activityIntegrate.setSleepPeriod(0);
        
        getController().addAction(activityIntegrate);
        integrator.setBox(box);
        
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
 
        //instantiate simulation
        NormalModeAnalysisDisplay3D sim = new NormalModeAnalysisDisplay3D(Space.getInstance(3));
        
        SimulationGraphic simGraphic = new SimulationGraphic(sim, Space.getInstance(3), sim.getController());
        simGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(50));
        simGraphic.makeAndDisplayFrame(APP_NAME);
        
	}
	
	public static class Applet extends javax.swing.JApplet{
		
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
	protected NormalModes nm;
	protected WaveVectorFactory waveVectorFactory;
	protected CoordinateDefinitionLeaf coordinateDefinition;
	protected ISpace space;
	protected double L;
	protected int n;
	
	protected static int numAtoms = 500;
	protected static double density = 1;
	protected static double temperature = 0.1;
	
}
