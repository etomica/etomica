package etomica.lattice;
import etomica.*;

/**
 * A 4-atom basis that makes an fcc crystal on a BravaisLattice
 * having a Cubic primitive.
 *
 * @author David Kofke
 */
 
 /* History
  * 09/22/02 (DAK) new
  */
 
public class BasisCubicFcc extends AtomFactoryHomo {
    
    /**
     * Makes a basis using a default that uses AtomFactoryMono
     * for making atom occupying each site.
     */
    public BasisCubicFcc(Simulation sim, PrimitiveCubic primitive) {
        this(sim, new AtomFactoryMono(sim), primitive);
    }
    /**
     * Makes a fcc 4-atom basis using the given factory to make the atoms.
     */
    public BasisCubicFcc(Simulation sim, AtomFactory factory, PrimitiveCubic primitive) {
        super(sim, factory, 4, BondInitializer.NULL, new ConfigurationCubicFcc(sim,primitive));
    }
    
    
    private static class ConfigurationCubicFcc extends Configuration {
        
        private ConfigurationCubicFcc(Simulation sim, PrimitiveCubic primitive) {
            super(sim);
            this.primitive = primitive;
        }
        
        private final Space3D.Vector[] positions = new Space3D.Vector[] {
            new Space3D.Vector(0.0, 0.0, 0.0),
            new Space3D.Vector(0.0, 0.5, 0.5),
            new Space3D.Vector(0.5, 0.5, 0.0),
            new Space3D.Vector(0.5, 0.0, 0.5)
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
    }//end ConfigurationCubicFcc
    
    
}//end of BasisCubicFcc