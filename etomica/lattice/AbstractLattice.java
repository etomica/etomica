package etomica.lattice;
import etomica.*;

public interface AbstractLattice extends java.io.Serializable {
    
    public int D();                      //dimension (1D, 2D, etc) of the lattice
    public AtomList siteList();          //list of all sites
    public Site site(int[] index);       //get a specific site
    public int[] dimensions();           //size of the lattice
    
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