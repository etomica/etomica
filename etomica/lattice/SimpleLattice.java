/*
 * History
 * Created on Dec 17, 2004 by kofke
 */
package etomica.lattice;

/**
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 *
 * @author David Kofke
 *
 */
public class SimpleLattice implements AbstractLattice {

    /**
     * 
     */
    public SimpleLattice(int D, SiteFactory siteFactory) {
        size = new int[D];
        size[0] = 1;
        dimensions = new int[D];
        this.siteFactory = siteFactory;
    }

    /* (non-Javadoc)
     * @see etomica.lattice.AbstractLattice#D()
     */
    public int D() {
        return dimensions.length;
    }

    /* (non-Javadoc)
     * @see etomica.lattice.AbstractLattice#siteList()
     */
    public Object[] sites() {
        return sites;
    }

    /* (non-Javadoc)
     * @see etomica.lattice.AbstractLattice#site(int[])
     */
    public Object site(int[] index) {
        int idx = 0;
        for(int i=0; i<dimensions.length; i++) {
            idx += index[i]*size[i];
        }
        return sites[idx];
    }

    /* (non-Javadoc)
     * @see etomica.lattice.AbstractLattice#getDimensions()
     */
    public int[] getDimensions() {
        return dimensions;
    }

    /* (non-Javadoc)
     * @see etomica.lattice.AbstractLattice#setDimensions(int[])
     */
    public void setDimensions(int[] dim) {
        if(dim.length != dimensions.length) throw new IllegalArgumentException("Incorrect dimension dimension");
        for(int i=0; i<dimensions.length; i++) {
            dimensions[i] = dim[i];
        }
        for(int i=0; i<size.length-1; i++) {
            size[i+1] = size[i]*dimensions[i];
        }
        sites = new Object[size[D()]*dimensions[D()]];
        siteFactory.makeSites(this, sites);
    }
    
    protected Object[] sites;
    protected final int[] dimensions;
    private final int[] size;
    protected SiteFactory siteFactory;
}
