/*
 * History
 * Created on Dec 17, 2004 by kofke
 */
package etomica.lattice;

import etomica.Space;


/**
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 *
 * @author David Kofke
 *
 */
public class CellLattice extends SimpleLattice {

    /**
     * @param D
     * @param siteFactory
     */
    public CellLattice(int D, SiteFactory siteFactory) {
        super(D, siteFactory);
        cellSize = new double[D];
        idx = new int[D];
    }
    
    public Object site(Space.Vector r) {
        for(int i=0; i<idx.length; i++) {
            idx[i] = (int)(r.x(i)/cellSize[i]);
//            double x = r.x(i)/cellSize[i];
//            idx[i] = (x < 0) ? (int)x - 1 : (int)x; //we want idx to be the floor of x
//            while(idx[i] >= dimensions[i]) idx[i] -= dimensions[i];
//            while(idx[i] < 0)              idx[i] += dimensions[i];
        }
        return site(idx);

    }

    public void setCellSize(double[] cellSize) {
        if(cellSize.length != this.cellSize.length) throw new IllegalArgumentException("Dimension of cell size array inconsistent with dimensions of cells");
        System.arraycopy(cellSize, 0, this.cellSize, 0, cellSize.length);
    }
    
    private final double[] cellSize;
    private final int[] idx;
}
