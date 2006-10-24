package etomica.eigenstuff;

import etomica.atom.AtomFactoryMono;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeSphere;
import etomica.config.ConfigurationLattice;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.models.hexane.PairIndexerMolecule;
import etomica.phase.Phase;
import etomica.potential.P2HardSphere;
import etomica.potential.Potential;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;

//Parts of this class are patterned on TestFcc.java
public class FccEigenGetter extends Simulation {

    //Constructor
    public FccEigenGetter(Space space, int numberOfAtoms, String filename){
        super(space, true, new PotentialMaster(space));
        dim = space.D();
        sysdim = numberOfAtoms * dim;
        
        LatticeCubicFcc lattice = new LatticeCubicFcc(1.0);
        ConfigurationLattice config = new ConfigurationLattice(lattice);
        // config.setRescalingToFitVolume(false);
        config.setRememberingIndices(true);

        defaults.makeLJDefaults();
        defaults.atomSize = 1.0;
        defaults.ignoreOverlap = true;

        SpeciesSpheresMono species = new SpeciesSpheresMono(this);

        phase = new Phase(this);
        phase.getAgent(species).setNMolecules(numberOfAtoms);

        config.initializeCoordinates(phase);
        
        //  *****Integrator stuff cut out.
        
        Potential potential = new P2HardSphere(space, defaults.atomSize, false);
        AtomTypeSphere sphereType = (AtomTypeSphere) ((AtomFactoryMono) species
                .moleculeFactory()).getType();
        potentialMaster.addPotential(potential, new AtomType[] { sphereType,
                sphereType });

        bdry = new BoundaryRectangularPeriodic(this);
        phase.setBoundary(bdry);

//      ***** more integrator stuff cut out
        
        config.initializeCoordinates(phase);

//      nan this section is a patch
        // first we find out the scaling used in
        // ConfigurationLattice/LatticeCubicFcc
        // then, we create a primitive fcc lattice, and scale it so we can use
        // it in pri.
        ConfigurationLattice.MyLattice myLattice = (ConfigurationLattice.MyLattice) config
                .getLatticeMemento();
        Vector scaling = myLattice.latticeScaling;
        scaling.TE(0.5 * Math.sqrt(2.0)); // we need this because the
                                            // fccPrimitive.setSize method uses
                                            // the unit vectors, not the actual
                                            // lattice vectors, and this
                                            // eliminates the radical 2 in the
                                            // unit vectors.
        Primitive fccPrimitive = lattice.getPrimitiveFcc();
        fccPrimitive.setSize(scaling.toArray());
        // nan phase.setDensity(1.04);
        // ***** Another integrator line cut out.
        
        pri = new PairIndexerMolecule(phase, fccPrimitive);
        
        // ***** output stuff cut out
        
        System.out.println(phase.getDensity());
        phase.setDensity(1.04);
        
        double[][][] fileData = ReaderIn.getFromFile(filename, dim);
        double[][] moreData = Reconstituter.reconstitute(pri, phase, fileData);
        evaler = new MyEigenvalueDecomposition(sysdim, moreData);
        double[] d = evaler.getRealEigenvalues();
        System.out.println("d  " + d);
        double[] e = evaler.getImagEigenvalues();
        System.out.println("e  " + e);
        double sum = 0;
        for (int i=0; i<d.length; i++) {
            sum += Math.log(2*Math.PI/d[i]);
//            System.out.println(d[i] + "   " + e[i] + "   " + sum );
        }
        System.out.println("f = "+sum);
        
    }
       
    
    
    
    public static void main (String[] args){
        int nA = 108;
        
        FccEigenGetter sim = new FccEigenGetter(Space3D.getInstance(), nA, "Happy.txt");
        
//        SimulationGraphic simG = new SimulationGraphic(sim);
//        simG.makeAndDisplayFrame();

        System.out.println("oik!");
    }
    
    MyEigenvalueDecomposition evaler;
    ReaderIn reader;
    Phase phase;
    BoundaryRectangularPeriodic bdry;
    PairIndexerMolecule pri;
    int dim;                            //the spatial dimension
    int sysdim;                         //the dimension of the matrix
}
