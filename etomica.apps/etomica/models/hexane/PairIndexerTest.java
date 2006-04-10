package etomica.models.hexane;

import nancyJunk.OutputFile;
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
import etomica.models.hexane.*;
import etomica.phase.Phase;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;
import etomica.util.Default;

public class PairIndexerTest extends Simulation {

    public PairIndexerTest(Space space, int numMolecules){
        super(space);
      
        Default def = new Default();
        def.makeLJDefaults();
        def.atomSize = 0.5;
        
        
        //****************  USE THIS FOR THE SIMPLE CUBIC TEST!!!!!
//        int chainLength = 1;
//        int numAtoms = numMolecules * chainLength;
//        Primitive prim = new PrimitiveCubic(space);
//        ConfigurationLattice config = new ConfigurationLattice(new BravaisLattice(prim));
//        config.setRememberingIndices(true);
//        Species species = new SpeciesSpheresMono(this);
//        species.setNMolecules(numMolecules);
//        def.boxSize = 6;
//        phase = new Phase(this);
//        bdry = new BoundaryRectangularPeriodic(space, def.boxSize);
        OutputFile printer = new OutputFile("Simple.txt");      
        
//        //****************  USE THIS FOR THE HEXANE TEST!!!!!
        int chainLength = 6;
        int numAtoms = numMolecules * chainLength;
        Primitive prim = new PrimitiveHexane(space);
        ConfigurationHexane config = new ConfigurationHexane(space);
        config.setRememberingIndices(true);
        Species species = new SpeciesHexane(this);
        species.setNMolecules(numMolecules);
        def.boxSize = 7.018;
        phase = new Phase(this);
        bdry = new BoundaryHexane(space);
        
//        OutputFile printer0 = new OutputFile("Printer0");
//        OutputFile printer1 = new OutputFile("Printer1");
//        OutputFile printer2 = new OutputFile("Printer2");
//        OutputFile printer3 = new OutputFile("Printer3");
//        OutputFile printer4 = new OutputFile("Printer4");
//        OutputFile printer5 = new OutputFile("Printer5");
        
        
        phase.setBoundary(bdry);
        config.initializeCoordinates(phase); 
        PairIndexer pi = new PairIndexer(phase, prim, chainLength);
        
        Atom atom0 = new AtomLeaf(phase.space());
        Atom atom1 = new AtomLeaf(phase.space());
        AtomPair ap = new AtomPair(atom0, atom1);
        
        AtomIterator inner = new AtomIteratorLeafAtoms(phase);
        AtomIterator outer = new AtomIteratorLeafAtoms(phase);
        ApiInnerFixed api = new ApiInnerFixed(outer, inner);
        api.reset();
        

        printer.println("Atom0:0 Atom0:1 Atom0:2 Atom0:# Atom1:0 Atom1:1 Atom1:2 Atom1:# Bin:");
        
//        printer0.print("Atom0:0 Atom0:1 Atom0:2 Atom0:# Atom1:0 Atom1:1 Atom1:2 Atom1:# Bin:");
//        printer1.print("Atom0:0 Atom0:1 Atom0:2 Atom0:# Atom1:0 Atom1:1 Atom1:2 Atom1:# Bin:");
//        printer2.print("Atom0:0 Atom0:1 Atom0:2 Atom0:# Atom1:0 Atom1:1 Atom1:2 Atom1:# Bin:");
//        printer3.print("Atom0:0 Atom0:1 Atom0:2 Atom0:# Atom1:0 Atom1:1 Atom1:2 Atom1:# Bin:");
//        printer4.print("Atom0:0 Atom0:1 Atom0:2 Atom0:# Atom1:0 Atom1:1 Atom1:2 Atom1:# Bin:");
//        printer5.print("Atom0:0 Atom0:1 Atom0:2 Atom0:# Atom1:0 Atom1:1 Atom1:2 Atom1:# Bin:");
        
        
        
        
        
        while(api.hasNext()){
            ap = api.nextPair();
            atom0 = ap.atom0;
            atom1 = ap.atom1;
            
            printer.println(pi.getIndex(atom0)[0] +" "+ pi.getIndex(atom0)[1] + " " +
                    pi.getIndex(atom0)[2] + " " + atom0.node.getIndex() + " " +
                    pi.getIndex(atom1)[0] + " " + pi.getIndex(atom1)[1] + " " +
                    pi.getIndex(atom1)[2] + " " + atom1.node.getIndex() + " " +
                    pi.getBin(ap));
   
        }
        
        printer.close();
       
//        printer0.close();
//        printer1.close();
//        printer2.close();
//        printer3.close();
//        printer4.close();
//        printer5.close();
        
    }
    
    
    /**
     * @param args
     */
    public static void main(String[] args) {
        int numMolecules = 216;
        
        PairIndexerTest pit = new PairIndexerTest(Space3D.getInstance(), numMolecules);
        
        System.out.println("Yes, we have no bananas!");
    
    }

    Phase phase;
    
    public Boundary bdry;
    
}
