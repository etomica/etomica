package etomica.lattice;


/**
 * Crystal built on a cubic lattice, such that its spatial size can be
 * characterized by a single value, the lattice constant.
 */

/*
 * History
 * Created on Jan 5, 2005 by kofke
 */
public interface LatticeCubic {
    
    public void setLatticeConstant(double latticeConstant);
    
    public double getLatticeConstant();

}
