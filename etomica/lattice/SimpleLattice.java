/*
 * History
 * Created on Dec 17, 2004 by kofke
 */
package etomica.lattice;

import java.awt.Color;
import java.awt.Graphics;

import javax.swing.JPanel;

import etomica.IteratorDirective;
import etomica.graphics.SimulationGraphic;

/**
 * Basic implementation of the AbstractLattice interface, providing construction
 * and access of sites for a lattice of arbitrary dimension.  Internally, sites are stored
 * in a 1-D array of objects, and are accessed by unrolling the index specification to
 * determine the storage-array index.
 */

//  D D D U U
//  D D D U U
//  D D X U U
//  D D U U U
//  D D U U U

public class SimpleLattice implements AbstractLattice {

    /**
     * Constructs a lattice of the given dimension (D) with sites
     * made from the given factory.
     */
    public SimpleLattice(int D, SiteFactory siteFactory) {
        this.D = D;
        jumpCount = new int[D];
        jumpCount[D-1] = 1;
        dimensions = new int[D];
        this.siteFactory = siteFactory;
    }

    /* (non-Javadoc)
     * @see etomica.lattice.AbstractLattice#D()
     */
    public final int D() {
        return D;
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
        return sites[arrayIndex(index)];
    }
    
    /**
     * Returns the index in the 1-d array for the site corresponding
     * to the given lattice index.
     */
    private final int arrayIndex(int [] index) {
        int idx = 0;
        for(int i=0; i<D; i++) {
            idx += index[i]*jumpCount[i];
        }
        return idx;
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
        if(dim.length != D) throw new IllegalArgumentException("Incorrect dimension dimension");
        System.arraycopy(dim, 0, dimensions, 0, D);
        for(int i=D-1; i>0; i--) {
            jumpCount[i-1] = jumpCount[i]*dimensions[i];
        }
        sites = new Object[jumpCount[0]*dimensions[0]];
        siteFactory.makeSites(this, sites);
    }
    
    /**
     * Returns a new instance of the neighbor iterator defined by the 
     * NeighborIterator inner class.
     */
    public NeighborIterator makeNeighborIterator() {
        return this.new NeighborIterator();
    }
    
    protected Object[] sites;
    protected final int[] dimensions;
//  jumpCount[i] gives the number of sites skipped when the i-th index is incremented by 1
    private final int[] jumpCount;
    private final int D;
    protected SiteFactory siteFactory;

    /**
     * An iterator that generates the neighboring sites of a given site.  The
     * neighbors are defined as those within a rectangular box centered on the 
     * site, with consideration of periodic boundaries.  Iterator can be configured
     * to yield the up- and/or down-neighbors of the site.  
     */
    public class NeighborIterator {

        /**
         * Constructs iterator that needs to be configured before use.  Must specify
         * the central site and the range before reset and iteration.
         */
        public NeighborIterator() {
            super();
            centralSite = new int[D];
            range = new int[D];
            cursor = Integer.MAX_VALUE;
        }
        
        /**
         * Identifies the site whose neighbors will be given on iteration.
         * Requires subsequent call to reset before iteration.
         * @param index the coordinate of the central site
         */
        public void setSite(int[] index) {
            if(index.length != D) throw new IllegalArgumentException("Incorrect length of array passed to setSite");
            System.arraycopy(index, 0, centralSite, 0, D);
            needNeighborUpdate = true;
            unset();
        }
        
        /**
         * Indicates whether the iterator has another site.
         */
        public final boolean hasNext() {
            return cursor < neighborCount;
        }
        
        /**
         * Resets the iterator to loop through its iterates again.  This
         * must be done after any call to setRange, setDirection, or setSite.
         */
        public void reset() {
            if(needNeighborUpdate) updateNeighborList();
            cursor = 0;
        }
        
        /**
         * Sets the range used to define the neighboring sites.  Each of the values
         * in newRange is used to determine the size of the neighbor box in the 
         * corresponding dimension: (edge length) = 2*value + 1, as the value describes
         * the neighbor distance in each direction (e.g., left and right) from the 
         * Non-negative values (including zero) are acceptable, to a maximum of
         * (lattice dimension - 1)/2, beyond which the neighbor distances extend beyond
         * the size of the lattice.  An IllegalArgumentException is thrown if any index
         * is outside its acceptable range.   
         */
        public void setRange(int[] newRange) {
            if(newRange.length != D) throw new IllegalArgumentException("Incorrect length of array passed to setRange");
            for(int i=0; i<newRange.length; i++) {
                if(newRange[i] < 0) 
                    throw new IllegalArgumentException("Neighbor range cannot be negative");
                if(2*newRange[i]+1 > dimensions[i]) 
                    throw new IllegalArgumentException("Neighbor range exceeds lattice site");
            }
            System.arraycopy(newRange, 0, range, 0, D);
            furthestNeighborDelta = arrayIndex(range);
            halfNeighborCount = 1;
            for(int i=0; i<D; i++) {
                halfNeighborCount *= (2*range[i]+1);
            }
            halfNeighborCount = (halfNeighborCount-1)/2;
            neighbors = new int[2*halfNeighborCount];
            needNeighborUpdate = true;
            unset();
        }

        /**
         * Method called on reset if any
         *
         */
        private void updateNeighborList() {
            neighborCount = (direction == null) ? 2*halfNeighborCount : halfNeighborCount;
            int centralSiteIndex = arrayIndex(centralSite);
            cursor = 0;
            if(doDown) gatherDownNeighbors(0, centralSiteIndex - furthestNeighborDelta);
            if(doUp) gatherUpNeighbors(0, centralSiteIndex+1);
            needNeighborUpdate = false;
        }
            
        private void gatherDownNeighbors(int d, int startIndex) {
            int centralIndex = centralSite[d];
            int iMin = centralIndex-range[d];
            int dim = dimensions[d];
            
            if(iMin < 0) {
                gatherNeighbors(d, -iMin, startIndex+dim*jumpCount[d]);
                gatherNeighbors(d, centralIndex, startIndex-iMin*jumpCount[d]);//note that iMin<0
            } else {
                gatherNeighbors(d, range[d], startIndex);
            }
            if(d < D-1) gatherDownNeighbors(d+1, startIndex+range[d]*jumpCount[d]);
        }

        private int gatherUpNeighbors(int d, int startIndex) {
            int centralIndex = centralSite[d];
            int iMax = centralIndex+range[d];
            int dim = dimensions[d];
            
            if(d < D-1) {
                startIndex = gatherUpNeighbors(d+1, startIndex);
                startIndex += jumpCount[d] - (range[d+1]+1)*jumpCount[d+1];
            }
            if(iMax >= dim) {
                gatherNeighbors(d, dim-centralIndex-1, startIndex);
                gatherNeighbors(d, iMax-dim+1, startIndex-(centralIndex+1)*jumpCount[d]);
            } else {
                gatherNeighbors(d, range[d], startIndex);
            }
            return startIndex;
        }
        
        private void gatherAllNeighbors(int d, int startIndex) {
            int centralIndex = centralSite[d];
            int iMin = centralIndex-range[d];
            int iMax = centralIndex+range[d];
            int dim = dimensions[d];
            
            if(iMin < 0) {
                gatherNeighbors(d, -iMin, startIndex+dim*jumpCount[d]);
                gatherNeighbors(d, iMax+1, startIndex-iMin*jumpCount[d]);//note that iMin<0
            } else if(iMax >= dim) {
                gatherNeighbors(d, dim-iMin, startIndex);
                gatherNeighbors(d, iMax-dim+1, startIndex-iMin*jumpCount[d]);
            } else {
                gatherNeighbors(d, 2*range[d]+1, startIndex);
            }
        }
        
        private void gatherNeighbors(int d, int nSteps, int startIndex) {
            if(d == D-1) {
                for(int i=0; i<nSteps; i++) {
                    neighbors[cursor++] = startIndex++;
                }
            } else {
                for(int i=0; i<nSteps; i++) {
                     gatherAllNeighbors(d+1,startIndex+i*jumpCount[d]);
                }
            }
        }
        
        /**
         * Puts iterator in a state in which hasNext() returns false.
         */
        public void unset() {
            cursor = neighborCount;
        }
        
        /**
         * Returns the next atom in the iteration sequence, or null
         * if hasNext is false.
         */
        public Object next() {
            return hasNext() ? sites[neighbors[cursor++]] : null;
        }
        
        /**
         * Returns the next iterate without advancing the iterator.
         */
        public Object peek() {
            return hasNext() ? sites[neighbors[cursor]] : null;
        }
        
        /**
         * The number of iterates returned by this iterator in its current state.
         */
        public int size() {
            return neighborCount;
        }

        /**
         * Gets the direction (up or down list) of the neighbors given by the iterator.
         * A null value indicates that all neighbors will be returned.
         */
       public IteratorDirective.Direction getDirection() {
            return direction;
        }
        /**
         * Sets the iterator to return only up-list neighbors or down-list neighbors.
         * A null value indicates that all neighbors are to be returned.
         */
        public void setDirection(IteratorDirective.Direction direction) {
            this.direction = direction;
            doUp = (direction != IteratorDirective.DOWN);//also handles case where
            doDown = (direction != IteratorDirective.UP);//direction is null
            needNeighborUpdate = true;
            unset();
        }
        
        private int[] range;
        private int[] neighbors;
        private final int[] centralSite;
        private int cursor;
        private int furthestNeighborDelta;
        private int neighborCount, halfNeighborCount;
        private boolean doUp, doDown;
        private boolean needNeighborUpdate = true;
        private IteratorDirective.Direction direction;
        private int siteIndex;
    }
  
    /**
     * Method to test the neighbor iterator.  Constructs a lattice 
     * with sites defined to hold a color.  Configures iterator, 
     * loops over neighbor sites and defines up/down neighbors to
     * have colors yellow and blue, respectively.  Then puts up a window
     * with the lattice drawn as an array of squares, with central site colored
     * red and with neighbor sites identified by their colors.  3-D lattice is shown
     * with a set of slices through the planes in the 3rd dimension.
     * @param args
     */
    public static void main(String[] args) {
        
        //define the site class
        class MySite {
            final int index;
            java.awt.Color color = java.awt.Color.white;
            MySite(int i) {this.index = i;}
        }
        //define a factory making the sites
        SiteFactory factory = new SiteFactory() {
            public void makeSites(AbstractLattice lattice, Object[] sites) {
                for(int i=0; i<sites.length; i++) {
                    sites[i] = new MySite(i);
                }
            }
        };
        //construct the lattice and iterator
        int dimension = 3;
        final SimpleLattice lattice = new SimpleLattice(dimension, factory);
        final SimpleLattice.NeighborIterator iterator = lattice.new NeighborIterator();
        //configure lattice
        //******change these values to perform different tests *********//
        //remember to change value of dimension in lattice constructor if using different dimensions here 
        lattice.setDimensions(new int[] {10,10,15});
        iterator.setSite(new int[] {8,5,12});
        iterator.setRange(new int[] {3,1,3});

        //define panel for display, so that it draws lattice with appropriately colored sites
        JPanel canvas = new JPanel() {
            public java.awt.Dimension getPreferredSize() {
                return new java.awt.Dimension(300,300);
            }
            public void paint(Graphics g) {
                int w = 8;  //size (pixels) of each square representing a site
                setSize(800,800);

                //generate down neighbors and color them yellow
                iterator.setDirection(IteratorDirective.DOWN);
                iterator.reset();
                while(iterator.hasNext()) {
                    MySite site = (MySite)iterator.next();
                    System.out.println(site.index);
                    site.color = java.awt.Color.YELLOW;
                }
                //generate up neighbors and color them blue
                iterator.setDirection(IteratorDirective.UP);
//                iterator.setDirection(null);
                iterator.reset();
                while(iterator.hasNext()) {
                    MySite site = (MySite)iterator.next();
                    System.out.println(site.index);
                    site.color = java.awt.Color.BLUE;
                }
                //color central site red
                ((MySite)lattice.site(iterator.centralSite)).color = Color.RED;
                int nx = 5; //number of planes to draw before moving down to another line of planes
                for(int k=0; k<lattice.dimensions[0]; k++) {//loop over planes
                    int ox = k % nx;       //set origin for drawing each lattice plane
                    int oy = (k - ox)/nx;
                    ox = 3 + ox*w*(lattice.dimensions[1]+2);
                    oy = 3 + oy*w*(lattice.dimensions[2]+2);
                    //draw lattice plane
                    for(int j=0; j<lattice.dimensions[2]; j++) {
                        for(int i=0; i<lattice.dimensions[1]; i++) {
                            g.setColor(((MySite)(lattice.site(new int[] {k,i,j}))).color);
                            g.fillRect(ox+i*w,oy+j*w,w,w);
                            g.setColor(Color.black);
                            g.drawRect(ox+i*w,oy+j*w,w,w);
                        }
                    }
                }
            }
        };
        //draw to screen
        javax.swing.JFrame f = new javax.swing.JFrame();
        f.setSize(700,500);
        f.getContentPane().add(canvas);
        f.pack();
        f.show();
        f.addWindowListener(SimulationGraphic.WINDOW_CLOSER);
    }//end of main
}
