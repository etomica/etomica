package etomica.lattice;
import etomica.AtomList;
import etomica.Atom;

public interface AbstractLattice extends java.io.Serializable {
    
    public int D();                      //dimension (1D, 2D, etc) of the lattice
    public AtomList siteList();          //list of all sites
    public Atom site(int[] index);       //get a specific site
    
    public interface Occupant {
        public Site site();
    }
    
}