package etomica.lattice;

/*
 * History of changes 09/18/02 (DAK) modified site method to return Atom instead
 * of Site.
 */
public interface AbstractLattice {

    public int D(); //dimension (1D, 2D, etc) of the lattice

    public Object[] sites(); //array of all sites

    public Object site(int[] index); //get a specific site

    public int[] getDimensions(); //size of the lattice

    public void setDimensions(int[] dim); //change lattice size to new value
}