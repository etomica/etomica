package etomica.lattice;

import java.awt.Color;
import java.awt.Graphics;

import javax.swing.JPanel;

import etomica.atom.iterator.IteratorDirective;
import etomica.graphics.SimulationGraphic;
import etomica.lattice.RectangularLattice.Iterator;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.Debug;


/**
 * Iterates over neighbors of a site, where neighbors are defined as
 * those sites lying in a rectangular region centered on the site.
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on May 26, 2005 by kofke
 */
public class RectangularLatticeNbrIteratorSquare extends
        RectangularLatticeNbrIterator {

    protected final int[] range;
    protected int furthestNeighborDelta;
    private final Vector period;
    protected int halfNeighborCount;

    /**
     * @param D
     */
    public RectangularLatticeNbrIteratorSquare(int D) {
        super(D);
        range = new int[D];
        Vector newPeriod = Space.makeVector(D);
        newPeriod.E(1.0);
        period = Space.makeVector(D);
        setPeriod(newPeriod);
    }

    /**
     * Method called on reset if any calls were previously made
     * to setLattice, setSite, setRange, or setDirection.
     */
    protected void updateNeighborList() {
        neighborCount = 0; //(direction == null) ? 2*halfNeighborCount : halfNeighborCount;
        int centralSiteIndex = lattice.arrayIndex(centralSite);
        cursor = 0;
        if(doDown) gatherDownNeighbors(0, centralSiteIndex - furthestNeighborDelta);
        if(doUp) gatherUpNeighbors(0, centralSiteIndex+1);
        needNeighborUpdate = false;
    }

    
    /* (non-Javadoc)
     * @see etomica.lattice.SiteIterator#setLattice(etomica.lattice.AbstractLattice)
     */
    public void setLattice(AbstractLattice lattice) {
        super.setLattice(lattice);
        for(int i=0; i<D; i++) {
            if(2*range[i]+1 > ((RectangularLattice)lattice).getSize()[i]) 
                 throw new IllegalArgumentException("Neighbor range exceeds lattice size");
         }

    }
    /**
     * Identifies and adds to neighbor list all downlist neighbors
     * in the sub-dimension d.  For example, if a 3D lattice and
     * d = 0, gathers all neighbors; d = 1, gathers all in a plane;
     * d = 2, gathers all in a line.  Neighbors are identified by their 
     * index in the sites array, and startIndex is the position in the array
     * where the first neighbor to be collected by the method is found.  
     */
    private void gatherDownNeighbors(int d, int startIndex) {
        int centralIndex = centralSite[d];
        int iMin = centralIndex-range[d];
        int dim = lattice.size[d];
        
        //this block gathers all neighbors in the dimension d,
        //up to but not including the row where the central site is located
        if(iMin < 0) {//need to implement periodic boundaries
            if (isPeriodic[d]) {
                gatherNeighbors(d, -iMin, startIndex+dim*lattice.jumpCount[d]);
            }
            gatherNeighbors(d, centralIndex, startIndex-iMin*lattice.jumpCount[d]);//note that iMin<0
        } else {//no concern for PBC; gather all neighbors in plane
            gatherNeighbors(d, range[d], startIndex);
        }
        //handle row containing central site
        if(d < D-1) gatherDownNeighbors(d+1, startIndex+range[d]*lattice.jumpCount[d]);
    }

    /**
     * Similar to gatherDownNeighbors, but does upList neighbors.
     */
    private int gatherUpNeighbors(int d, int startIndex) {
        int centralIndex = centralSite[d];
        int iMax = centralIndex+range[d];
        int dim = lattice.size[d];
        
        //handle central-site row
        if(d < D-1) {
            startIndex = gatherUpNeighbors(d+1, startIndex);
            startIndex += lattice.jumpCount[d] - (range[d+1]+1)*lattice.jumpCount[d+1];
        }
        //advance through other rows
        if(iMax >= dim) {
            gatherNeighbors(d, dim-centralIndex-1, startIndex);
            if (isPeriodic[d]) {
                gatherNeighbors(d, iMax-dim+1, startIndex-(centralIndex+1)*lattice.jumpCount[d]);
            }
        } else {
            gatherNeighbors(d, range[d], startIndex);
        }
        return startIndex;
    }

    /**
     * @param d specifies row, plane, block, etc. of neighbors gathered
     * @param nSteps number of steps to take in row
     * @param startIndex array of site index for first neighbor
     */
    private void gatherNeighbors(int d, int nSteps, int startIndex) {
        if(d == D-1) {//end of recursion -- here's where the actual gathering is done
            for(int i=0; i<nSteps; i++) {
                neighborCount++;
                neighbors[cursor++] = startIndex++;
            }
        } else {//step from one row to the next and gather neighbors in each row
            int d1 = d+1;
            int centralIndex = centralSite[d1];
            int iMin = centralIndex-range[d1];
            int iMax = centralIndex+range[d1];
            int dim = lattice.size[d1];
            if(iMin < 0) {//need to consider PBC below
                for(int i=0; i<nSteps; i++) {
                    int startIndex1 = startIndex+i*lattice.jumpCount[d];
                    if (isPeriodic[d1]) {
                        gatherNeighbors(d1, -iMin, startIndex1+dim*lattice.jumpCount[d1]);
                    }
                    gatherNeighbors(d1, iMax+1, startIndex1-iMin*lattice.jumpCount[d1]);//note that iMin<0
                }
            } else if(iMax >= dim) {//need to consider PBC above
                for(int i=0; i<nSteps; i++) {
                    int startIndex1 = startIndex+i*lattice.jumpCount[d];
                    gatherNeighbors(d1, dim-iMin, startIndex1);
                    if (isPeriodic[d1]) {
                        gatherNeighbors(d1, iMax-dim+1, startIndex1-iMin*lattice.jumpCount[d1]);
                    }
                }
            } else {//no need to consider PBC
                for(int i=0; i<nSteps; i++) {
                    gatherNeighbors(d1, 2*range[d1]+1, startIndex+i*lattice.jumpCount[d]);
                }
            }
        }
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
     * A copy of the newRange array is made, so it may be re-used elsewhere after
     * calling this method.
     */
    public void setRange(int[] newRange) {
        if(Debug.ON) {
            if(newRange.length != D) throw new IllegalArgumentException("Incorrect length of array passed to setRange");
            for(int i=0; i<D; i++) {
                if(newRange[i] < 0) 
                    throw new IllegalArgumentException("Neighbor range cannot be negative");
                if(2*newRange[i]+1 > lattice.size[i]) 
                    throw new IllegalArgumentException("Neighbor range exceeds lattice site");
            }
        }
        halfNeighborCount = 1;
        for(int i=0; i<D; i++) {
            range[i] = newRange[i];
            halfNeighborCount *= (2*range[i]+1);
        }
        halfNeighborCount = (halfNeighborCount-1)/2;
        neighbors = new int[2*halfNeighborCount];
        furthestNeighborDelta = lattice.arrayIndex(range);
        needNeighborUpdate = true;
        unset();
    }

    /**
     * Returns (a copy of) the array specifying the neighbor range.
     */
    public int[] getRange() {
        return (int[])range.clone();
    }

    /**
     * Sets the magnitude of the elements of the pbc array.  Normally
     * this would be the dimensions of the phase boundary.  This quantity
     * has no effect on the selection of the neighbor cells; it affects only
     * the vector returned by the currentPbc method.
     * @param newPeriod values of new period are copied to internal vector
     */
    public void setPeriod(Vector newPeriod) {
        if(period.equals(newPeriod)) return;
        
        period.E(newPeriod);
    
        int[] idx = new int[D];
        for (int i=0; i<D; i++) {
            idx[i] = -1;
        }
        double[] vectorElements = new double[D];
        // increments the last dimension first, starting at {-1,-1,-1}
        cursorJump[D-1] = 1;
        for (int i=D-2; i>-1; i--) {
            cursorJump[i] = cursorJump[i+1]*3;
        }
        while (idx[0] < 2) {
            for(int j=0; j<D; j++) {
                vectorElements[j] = idx[j]*period.x(j);
            }
            int i=D-1;
            idx[i]++;
            while (idx[i] > 1 && i > 0) {
                idx[i] = -1;
                idx[--i]++;
            }
        }
    
    }

//  private void gatherAllNeighbors(int d, int startIndex) {
//  int centralIndex = centralSite[d];
//  int iMin = centralIndex-range[d];
//  int iMax = centralIndex+range[d];
//  int dim = dimensions[d];
//  
//  if(iMin < 0) {//need to consider PBC
//      gatherNeighbors(d, -iMin, startIndex+dim*jumpCount[d]);
//      gatherNeighbors(d, iMax+1, startIndex-iMin*jumpCount[d]);//note that iMin<0
//  } else if(iMax >= dim) {//need to consider PBC
//      gatherNeighbors(d, dim-iMin, startIndex);
//      gatherNeighbors(d, iMax-dim+1, startIndex-iMin*jumpCount[d]);
//  } else {//no need to consider PBC
//      gatherNeighbors(d, 2*range[d]+1, startIndex);
//  }
//}
//
///**
//* @param d specifies row, plane, block, etc. of neighbors gathered
//* @param nSteps number of steps to take in row
//* @param startIndex array of site index for first neighbor
//*/
//private void gatherNeighbors(int d, int nSteps, int startIndex) {
//  if(d == D-1) {//end of recursion -- here's where the actual gathering is done
//      for(int i=0; i<nSteps; i++) {
//          neighbors[cursor++] = startIndex++;
//      }
//  } else {//step from one row to the next and gather neighbors in each row
//      for(int i=0; i<nSteps; i++) {
//           gatherAllNeighbors(d+1,startIndex+i*jumpCount[d]);
//      }
//  }
//}

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
            final int[] coord;
            java.awt.Color color = java.awt.Color.white;
            MySite(int i, int[] idx) {this.index = i;this.coord=(int[])idx.clone();}
        }
        //define a factory making the sites
        SiteFactory factory = new SiteFactory() {
            public Object makeSite(AbstractLattice lattice, int[] coord) {
                return new MySite(((RectangularLattice)lattice).arrayIndex(coord),coord);
            }
//            public void makeSites(AbstractLattice lattice, Object[] sites) {
//                for(int i=0; i<sites.length; i++) {
//                    sites[i] = new MySite(i);
//                }
//            }
        };
        //construct the lattice and iterator
        int dimension = 3;
        final RectangularLattice lattice = new RectangularLattice(dimension, factory);
        final RectangularLattice.Iterator siteIterator = new Iterator(dimension);
        siteIterator.setLattice(lattice);
        //configure lattice
        //******change these values to perform different tests *********//
        //remember to change value of dimension in lattice constructor if using different dimensions here 
        lattice.setSize(new int[] {12,11,15});

        final RectangularLatticeNbrIteratorSquare iterator = new RectangularLatticeNbrIteratorSquare(dimension);
        iterator.setPeriodicity(new boolean[]{true,true,false});
        iterator.setLattice(lattice);
        iterator.setSite(new int[] {1,1,14});
        iterator.setRange(new int[] {2,3,5});
        //define panel for display, so that it draws lattice with appropriately colored sites
        JPanel canvas = new JPanel() {
            public java.awt.Dimension getPreferredSize() {
                return new java.awt.Dimension(600,500);
            }
            public void paint(Graphics g) {
                int w = 8;  //size (pixels) of each square representing a site
                setSize(800,800);

                //generate down neighbors and color them green
                iterator.setDirection(IteratorDirective.Direction.DOWN);
                iterator.reset();
                while(iterator.hasNext()) {
                    MySite site = (MySite)iterator.next();
//                    System.out.println(site.index);
                    site.color = java.awt.Color.GREEN.darker();
                }
                //generate up neighbors and color them blue
                iterator.setDirection(IteratorDirective.Direction.UP);
//                iterator.setDirection(null);
                iterator.reset();
                while(iterator.hasNext()) {
                    MySite site = (MySite)iterator.next();
//                    System.out.println(site.index);
                    site.color = java.awt.Color.BLUE.darker();
                }
//                siteIterator.reset();
//                while(siteIterator.hasNext()) {
//                    MySite site = (MySite)siteIterator.peek();
//                    int[] coord = (int[])siteIterator.nextIndex();
//                    System.out.print(site.index + "[");
//                    for(int i=0; i<lattice.D(); i++) System.out.print(site.coord[i]+" ");
//                    System.out.print(" [");
//                    for(int i=0; i<lattice.D(); i++) System.out.print(coord[i]+" ");
//                    System.out.print(" [");
//                    int[] latticeIndex = lattice.latticeIndex(site.index);
//                    for(int i=0; i<lattice.D(); i++) System.out.print(latticeIndex[i]+" ");
//                    System.out.println();
//                    site.color = java.awt.Color.BLUE;
//                }
                
                //color central site red
//                ((MySite)lattice.site(iterator.centralSite)).color = Color.RED;
                int nx = 5; //number of planes to draw before moving down to another line of planes
                for(int k=0; k<lattice.size[0]; k++) {//loop over planes
                    int ox = k % nx;       //set origin for drawing each lattice plane
                    int oy = (k - ox)/nx;
                    ox = 3 + ox*w*(lattice.size[1]+2);
                    oy = 3 + oy*w*(lattice.size[2]+2);
                    //draw lattice plane
                    for(int j=0; j<lattice.size[2]; j++) {
                        for(int i=0; i<lattice.size[1]; i++) {
                            g.setColor(((MySite)(lattice.site(new int[] {k,i,j}))).color);
                            g.fillRect(ox+i*w,oy+j*w,w,w);
                            g.setColor(Color.black);
                            g.drawRect(ox+i*w,oy+j*w,w,w);
                        }
                    }
                }
            }//end of paint
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
