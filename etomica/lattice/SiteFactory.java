package etomica.lattice;

/**
 * Interface for classes that can construct Site objects.
 *
 * @author David Kofke
 */
public interface SiteFactory extends java.io.Serializable {
    
   public Site makeSite(AbstractLattice parent, SiteIterator.Neighbor iterator, AbstractLattice.Coordinate coord);

}