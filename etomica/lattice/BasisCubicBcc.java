package etomica.lattice;
import etomica.*;

/**
 * A 2-atom basis that makes a bcc crystal on a BravaisLattice
 * having a Cubic primitive.
 *
 * @author David Kofke
 */
 
 /* History
  * 09/22/02 (DAK) new
  */
 
public class BasisCubicBcc extends AtomFactoryHomo {
    
    /**
     * Makes a basis using a default that uses AtomFactoryMono
     * for making atom occupying each site.
     */
    public BasisCubicBcc(Simulation sim, PrimitiveCubic primitive) {
        this(sim, new AtomFactoryMono(sim), primitive);
    }
    /**
     * Makes a bcc 2-atom basis using the given factory to make the atoms.
     */
    public BasisCubicBcc(Simulation sim, AtomFactory factory, PrimitiveCubic primitive) {
        super(sim, factory, 2, BondInitializer.NULL, new ConfigurationCubicBcc(sim,primitive));
    }
    
    
    private static class ConfigurationCubicBcc extends Configuration {
        
        private ConfigurationCubicBcc(Simulation sim, PrimitiveCubic primitive) {
            super(sim);
            this.primitive = primitive;
        }
        
        private final Space3D.Vector[] positions = new Space3D.Vector[] {
            new Space3D.Vector(0.0, 0.0, 0.0),
            new Space3D.Vector(0.5, 0.5, 0.5)
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
                iterator.next().coord.translateTo(r);
            }
        }
    }//end ConfigurationCubicBcc
    
    
}//end of BasisCubicBcc