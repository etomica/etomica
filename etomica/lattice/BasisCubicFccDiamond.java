package etomica.lattice;
import etomica.*;

/**
 * A 2-atom basis that makes a diamond crystal using a BravaisLattice
 * having a Cubic primitive with an fcc basis.  Each of the 4 fcc sites
 * is populated by a 2-atom basis arranged to yield the diamond structure.
 *
 * @author David Kofke
 */
 
 /* History
  * 09/26/02 (DAK) new
  * 01/19/04 (DAK) revised to extend Basis instead of AtomFactory
  */
 
public class BasisCubicFccDiamond extends Basis {
    
    /**
     * Makes a basis using a default that uses AtomFactoryMono
     * for making atom occupying each site.
     * @param space instance of governing space class
     * @param primitive Primitive of the cubic lattice housing this basis.
     * Needed to ensure that separation of basis atoms is consistent with
     * spacing of atoms on lattice.     */
    public BasisCubicFccDiamond(Space space, PrimitiveCubic primitive) {
		super(space, 2, new Configuration(space,primitive));
    }
    /**
     * Makes a diamond-on-fcc 2-atom basis using the given factory to make the atoms.
	 * @param space 
	 * @param size Number of atoms (n) in basis
	 * @param configuration Specifies arrangement of atoms
	 * @param atomFactory Used to constructe each site in the basis
     */
//    public BasisCubicFccDiamond(Space space, AtomSequencer.Factory seqFactory, 
//                                AtomFactory factory, PrimitiveCubic primitive) {
//        super(space, seqFactory, factory, 2, BondInitializer.NULL, new Configuration(space,primitive));
	public BasisCubicFccDiamond(Space space, AtomFactory factory, PrimitiveCubic primitive) {
		super(space, 2, new Configuration(space,primitive), factory);
    }
    
    
    private static class Configuration extends etomica.Configuration {
        
        private Configuration(Space space, PrimitiveCubic primitive) {
            super(space);
            this.primitive = primitive;
        }
        
        private final Space3D.Vector[] positions = new Space3D.Vector[] {
            new Space3D.Vector(0.0, 0.0, 0.0),
            new Space3D.Vector(0.25, 0.25, 0.25),
        };
        private final Space3D.Vector r = new Space3D.Vector();
        private PrimitiveCubic primitive;
        
            
        public void initializePositions(AtomIterator[] iterators){
            if(iterators == null || iterators.length == 0) return;
            AtomIterator iterator;
            if(iterators.length == 1) iterator = iterators[0];
            else iterator = new AtomIteratorCompound(iterators);//lump 'em all together
            iterator.reset();
            double latticeConstant = primitive.getSize();
            int i = 0;
            while(iterator.hasNext()) {
                r.Ea1Tv1(latticeConstant,positions[i++]);
                Atom a = iterator.nextAtom();
                try {//may get null pointer exception when beginning simulation
                    a.creator().getConfiguration().initializePositions(a);
                } catch(NullPointerException e) {}
                a.coord.translateTo(r);
            }
        }
    }//end Configuration
    
    
}//end of BasisCubicFcc