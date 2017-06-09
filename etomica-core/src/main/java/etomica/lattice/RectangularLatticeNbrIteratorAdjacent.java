/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.lattice;

import etomica.graphics.SimulationGraphic;
import etomica.lattice.RectangularLattice.Iterator;
import etomica.potential.IteratorDirective;

import javax.swing.*;
import java.awt.*;


/**
 * Iterates over the neighbors of a site, where neighbors are those
 * sites immediately adjacent to it in each direction; this does not
 * include sites diagonal to the site.  In D dimensions, there are 2D
 * such neighbors.  The neighbor definition may or may not be defined
 * as obeying periodicity, as desired.
 *
 * @author David Kofke
 *
 */
public class RectangularLatticeNbrIteratorAdjacent extends
        RectangularLatticeNbrIterator {


    /**
     * @param D
     */
    public RectangularLatticeNbrIteratorAdjacent(int D) {
        super(D);
        neighbors = new int[2*D];
    }

    /* (non-Javadoc)
     * @see etomica.lattice.RectangularLatticeNbrIterator#updateNeighborList()
     */
    protected void updateNeighborList() {
        neighborCount = 0; //(direction == null) ? 2*halfNeighborCount : halfNeighborCount;
        int centralSiteIndex = lattice.arrayIndex(centralSite);
        cursor = 0;
        if(doDown) {
            for(int d=0; d<D; d++) {
                if(centralSite[d] == 0) {
                    if(isPeriodic[d]) {
                        int idx = centralSiteIndex - lattice.jumpCount[d];
                        if(d > 0) neighbors[cursor++] = idx + lattice.jumpCount[d-1];
                        else neighbors[cursor++] = idx + lattice.sites.length;
                        neighborCount++;
                    }
                }
                else {
                    neighbors[cursor++] = centralSiteIndex - lattice.jumpCount[d];
                    neighborCount++;
                }
            }
        }
        if(doUp) {
            for(int d=0; d<D; d++) {
                if(centralSite[d] == lattice.size[d]-1) {
                    if(isPeriodic[d]) {
                        int idx = centralSiteIndex + lattice.jumpCount[d];
                        if(d > 0) neighbors[cursor++] = idx - lattice.jumpCount[d-1];
                        else neighbors[cursor++] = idx - lattice.sites.length;
                        neighborCount++;
                    }
                }
                else {
                    neighbors[cursor++] = centralSiteIndex + lattice.jumpCount[d];
                    neighborCount++;
                }
            }
        }

        needNeighborUpdate = false;
    }
    
    private static final long serialVersionUID = 1L;

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
            java.awt.Color color = java.awt.Color.white;
        }
        //define a factory making the sites
        SiteFactory factory = new SiteFactory() {
            public Object makeSite(AbstractLattice lattice, int[] coord) {
                return new MySite();
            }
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

        final RectangularLatticeNbrIteratorAdjacent iterator = new RectangularLatticeNbrIteratorAdjacent(dimension);
        iterator.setPeriodicity(new boolean[]{true,true,true});
        iterator.setLattice(lattice);
        iterator.setSite(new int[] {5,2,7});
        //define panel for display, sothat it draws lattice with appropriately colored sites
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
                
                //color central site red
                ((MySite)lattice.site(iterator.centralSite)).color = Color.RED;
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
    }
}
