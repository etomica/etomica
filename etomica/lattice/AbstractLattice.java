package etomica.lattice;
import etomica.*;

/* History of changes
 * 09/18/02 (DAK) modified site method to return Atom instead of Site.
 */
public interface AbstractLattice extends java.io.Serializable {
    
    public int D();                      //dimension (1D, 2D, etc) of the lattice
    public AtomList siteList();          //list of all sites
    public Atom site(int[] index);       //get a specific site
    public int[] getDimensions();           //size of the lattice
    public void setDimensions(int[] dim); //change lattice size to new value
    public SimulationEventManager eventManager(); //manages listeners and event firing
    
    public interface Occupant {
        public Site site();
    }
    
    
    //consider this
    
/*    public interface LatticeNode {
        public boolean isSite(); //distinguish a lattice site from a node in the tree
        public AbstractLattice parentLattice();
    }
    public class SiteNode extends AtomTreeNodeLeaf implements LatticeNode {}
    public class GroupNode extends AtomTreeNodeGroup implements LatticeNode {}
    */
}