package etomica.normalmode;

import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomList;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.api.IVectorMutable;
import etomica.atom.iterator.IteratorDirective;
import etomica.box.Box;
import etomica.integrator.IntegratorMC;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.potential.P2SoftSphere;
import etomica.potential.P2SoftSphericalTruncatedShifted;
import etomica.potential.Potential2SoftSpherical;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;

/**
 * 
 *	
 * 
 * @author Tai Boon Tan
 */
public class EinsteinCrystalExpansionAnalysis extends Simulation {

    public EinsteinCrystalExpansionAnalysis(Space _space, int numAtoms, double density, int exponent) {
        super(_space, true);


        potentialMaster = new PotentialMasterMonatomic(this);

        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        getSpeciesManager().addSpecies(species);

        box = new Box(space);
        addBox(box);
        box.setNMolecules(species, numAtoms);

        
        double L = Math.pow(4.0/density, 1.0/3.0);
        primitive = new PrimitiveCubic(space, L);
        int n = (int)Math.round(Math.pow(numAtoms/4, 1.0/3.0));
        nCells = new int[]{n,n,n};
        boundary = new BoundaryRectangularPeriodic(space, n * L);
        box.setBoundary(boundary);
        basis = new BasisCubicFcc();
        
        coordinateDefinition = new CoordinateDefinitionLeaf(this, box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(nCells);
        
        IVectorMutable[] initialLatticePos = _space.makeVectorArray(numAtoms);
        IAtomList atoms = box.getLeafList();
        
        for (int i=0; i<numAtoms; i++){
        	initialLatticePos[i].E(coordinateDefinition.getLatticePosition(atoms.getAtom(i)));
        }
      
      /*
       *   Reassigning the density
       */
        
        energySum = new PCEnergySumEinsteinCrystalExpansion();
        energySum.setInitialLatticePosition(initialLatticePos);
        energySum.setBox(box);
        IteratorDirective id = new IteratorDirective();
        
        Potential2SoftSpherical potential = new P2SoftSphere(space);
        P2SoftSphericalTruncatedShifted pTruncated = new P2SoftSphericalTruncatedShifted(space, potential, 1);
	    
        IAtomType sphereType = species.getLeafType();
        potentialMaster.addPotential(pTruncated, new IAtomType[] {sphereType, sphereType});
        
        int looping = 10000;
        double interval = 0.0001;
        double[] value = new double[looping];
        double gamma =1.0;
        
        for (int i=0; i<looping; i++){
        
        	
        	
	        double reducedDensity = density / (gamma* gamma* gamma);
	        L = Math.pow(4.0/reducedDensity, 1.0/3.0);
	        primitive = new PrimitiveCubic(space, L);
	        
	        boundary.setBoxSize(_space.makeVector(new double[]{n*L, n*L, n*L}));
	        
	        coordinateDefinition = new CoordinateDefinitionLeaf(this, box, primitive, basis, space);
	        coordinateDefinition.initializeCoordinates(nCells);
	        /*
	         * 
	         */
	        
	        double truncationRadius = boundary.getBoxSize().getX(0) * 0.495;
	        pTruncated.setTruncationRadius(truncationRadius);

	        energySum.zeroSum();
	        potentialMaster.calculate(box, id, energySum);
	        
	        value[i] = energySum.getSum()/numAtoms;
	        
	        System.out.println("The value["+i+"] is: "+ value[i] + ", gamma: "+gamma);
	        gamma += interval;
	        
        } //end looping over gamma
        
        System.out.println("The sum value is: "+trapezoidalIntegration(value, interval, looping));
        
    }

    
    public double trapezoidalIntegration(double[] f, double interval, int looping){
    	
    	if (looping ==1){
    		return f[0];
    	}
    	
    	double value = 0;
    	for (int i=1; i<looping; i++){
    		value += 0.5*(f[i]+f[i-1])*interval;
    	}
    	return value;
    }
    /**
     * @param args
     */
    public static void main(String[] args) {

        // defaults
        int D = 3;
        int nA = 32;
        double density = 1.256;
        int exponent = 12 ;


        // construct simulation
        EinsteinCrystalExpansionAnalysis sim = new EinsteinCrystalExpansionAnalysis(Space.getInstance(D), nA, density, exponent);


      
  }

    private static final long serialVersionUID = 1L;
    public IntegratorMC integrator;
    public ActivityIntegrate activityIntegrate;
    public IBox box;
    public Boundary boundary;
    public Primitive primitive;
    public Basis basis;
    public int[] nCells;
    public CoordinateDefinition coordinateDefinition;
    public PotentialMasterMonatomic potentialMaster;
    public PCEnergySumEinsteinCrystalExpansion energySum;
}