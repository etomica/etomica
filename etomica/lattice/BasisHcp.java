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
  * 01/19/04 (DAK) revised to extend Basis instead of AtomFactory
  */
 
public class BasisHcp extends Basis {
    
    /**
     * Makes a basis using a default that uses AtomFactoryMono
     * for making atom occupying each site.
     * @param space instance of governing space class
     * @param primitive Primitive of the cubic lattice housing this basis.
     * Needed to ensure that separation of basis atoms is consistent with
     * spacing of atoms on lattice.
     */
    public BasisHcp(Space space, PrimitiveHexagonal primitive) {
		super(space, 2, new Configuration(space,primitive));
    }
    /**
     * Makes a hcp 2-atom basis using the given factory to make the atoms.
     * @param space instance of governing space class
     * @param factory AtomFactory used to make the atoms (molecules) the form
     * the basis
     * @param primitive Primitive of the cubic lattice housing this basis.
     * Needed to ensure that separation of basis atoms is consistent with
     * spacing of atoms on lattice.
      */
    public BasisHcp(Space space, AtomFactory factory, PrimitiveHexagonal primitive) {
		super(space, 2, new Configuration(space,primitive), factory);
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
            Atom atom = iterator.nextAtom(); 
            try {//may get null pointer exception when beginning simulation
                atom.creator().getConfiguration().initializePositions(atom);
            } catch(NullPointerException e) {}
            atom.coord.translateTo(r);
            
            //second atom
            for(int i=0; i<a.length; i++) {
                r.PEa1Tv1(factors[i], a[i]);
            }
            atom = iterator.nextAtom();
            try {//may get null pointer exception when beginning simulation
                atom.creator().getConfiguration().initializePositions(atom);
            } catch(NullPointerException e) {}
            atom.coord.translateTo(r);
        }
    }//end Configuration
    
    
}//end of BasisHcp