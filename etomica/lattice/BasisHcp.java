package etomica.lattice;
import etomica.*;

/**
 * A 2-atom basis that makes a hcp crystal on a BravaisLattice
 * having a Hexagonal primitive.
 *
 * @author David Kofke
 */
 
 /* History
  * 09/22/02 (DAK) new
  */
 
public class BasisHcp extends AtomFactoryHomo {
    
    /**
     * Makes a basis using a default that uses AtomFactoryMono
     * for making atom occupying each site.
     */
    public BasisHcp(Space space, PrimitiveHexagonal primitive) {
        this(space, new AtomFactoryMono(space, AtomSequencerSimple.FACTORY), primitive);
    }
    /**
     * Makes a hcp 2-atom basis using the given factory to make the atoms.
     */
    public BasisHcp(Space space, AtomFactory factory, PrimitiveHexagonal primitive) {
        super(space, AtomSequencerSimple.FACTORY, factory, 2, BondInitializer.NULL, new Configuration(space,primitive));
    }
    
    
    private static class Configuration extends etomica.Configuration {
        
        private Configuration(Space space, PrimitiveHexagonal primitive) {
            super(space);
            this.primitive = primitive;
        }
        
        private final double[] factors = new double[] {1./3., 1./3., 0.5};
        private final Space3D.Vector r = new Space3D.Vector();
        private PrimitiveHexagonal primitive;
        
            
        public void initializePositions(AtomIterator[] iterators){
            if(iterators == null || iterators.length == 0) return;
            AtomIterator iterator;
            if(iterators.length == 1) iterator = iterators[0];
            else iterator = new AtomIteratorCompound(iterators);//lump 'em all together
            Space.Vector[] a = primitive.vectors();
            iterator.reset();
            //first atom
            r.E(0.0);
            Atom atom = iterator.next(); 
            try {//may get null pointer exception when beginning simulation
                atom.creator().getConfiguration().initializePositions(atom);
            } catch(NullPointerException e) {}
            atom.coord.translateTo(r);
            
            //second atom
            for(int i=0; i<a.length; i++) {
                r.PEa1Tv1(factors[i], a[i]);
            }
            atom = iterator.next();
            try {//may get null pointer exception when beginning simulation
                atom.creator().getConfiguration().initializePositions(atom);
            } catch(NullPointerException e) {}
            atom.coord.translateTo(r);
        }
    }//end Configuration
    
    
}//end of BasisHcp