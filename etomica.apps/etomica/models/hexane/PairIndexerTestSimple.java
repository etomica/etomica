package etomica.models.hexane;

import etomica.atom.Atom;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomPair;
import etomica.atom.iterator.ApiInnerFixed;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.config.ConfigurationLattice;
import etomica.lattice.BravaisLattice;
import etomica.lattice.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.phase.Phase;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;
import etomica.util.Default;

public class PairIndexerTestSimple extends Simulation {

    public PairIndexerTestSimple(Space space, int numMolecules){
        super(space);
      
        Default def = new Default();
        def.makeLJDefaults();
        def.atomSize = 0.5;
        
        int chainLength = 1;
        int numAtoms = numMolecules * chainLength;
        Primitive prim = new PrimitiveCubic(space);
        ConfigurationLattice config = new ConfigurationLattice(new BravaisLattice(prim));
        config.setRememberingIndices(true);
        Species species = new SpeciesSpheresMono(this);
        species.setNMolecules(numMolecules);
        def.boxSize = 6;
        phase = new Phase(this);
        bdry = new BoundaryRectangularPeriodic(space, def.boxSize);
        OutputFile printer = new OutputFile("Simple.txt");      
        
        phase.setBoundary(bdry);
        config.initializeCoordinates(phase); 
        //nan this will need to be changed
        PairIndexerMolecule pi = new PairIndexerMolecule(phase, prim);
        
        Atom atom0 = new AtomLeaf(phase.space());
        Atom atom1 = new AtomLeaf(phase.space());
        AtomPair ap = new AtomPair(atom0, atom1);
        
        AtomIterator inner = new AtomIteratorLeafAtoms(phase);
        AtomIterator outer = new AtomIteratorLeafAtoms(phase);
        ApiInnerFixed api = new ApiInnerFixed(outer, inner);
        api.reset();
        

        printer.println("Atom0:0 Atom0:1 Atom0:2 Atom0:# Atom1:0 Atom1:1 Atom1:2 Atom1:# Bin:");
        
        while(api.hasNext()){
            ap = api.nextPair();
            atom0 = ap.atom0;
            atom1 = ap.atom1;
        
            printer.println(pi.getIndex(atom0)[0] +" "+ pi.getIndex(atom0)[1] + " " +
                    pi.getIndex(atom0)[2] + " " + atom0.node.getIndex() + " " +
                    pi.getIndex(atom1)[0] + " " + pi.getIndex(atom1)[1] + " " +
                    pi.getIndex(atom1)[2] + " " + atom1.node.getIndex() + " " +
                    pi.getBin(ap) +" " + atom0.getGlobalIndex());
   
        }
        
        printer.close();
          
    }
    
    
    /**
     * @param args
     */
    public static void main(String[] args) {
        int numMolecules = 216;
        
        PairIndexerTestSimple pit = new PairIndexerTestSimple(Space3D.getInstance(), numMolecules);
        
        System.out.println("Yes, we have no bananas!");
    
    }

    Phase phase;
    
    public Boundary bdry;
    
}
