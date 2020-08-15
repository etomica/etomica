/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;

import java.awt.Color;
import java.awt.Graphics;

import etomica.action.controller.Controller;
import etomica.atom.IAtom;
import etomica.lattice.RectangularLattice;
import etomica.nbr.site.AtomSite;
import etomica.nbr.site.NeighborSiteManager;
import etomica.space.Space;


public class DisplayBoxSpin2D extends DisplayBoxCanvas2D {

    /**
     * @param _box
     */
    public DisplayBoxSpin2D(DisplayBox _box, NeighborSiteManager neighborSiteManager, Space space, Controller controller) {
        super(_box, space, controller);
        latticeIndex = new int[displayBox.getBox().getBoundary().getBoxSize().getD()];
        spinWidth = 5;
        this.neighborSiteManager = neighborSiteManager;
    }
    
    protected void drawAtom(Graphics g, int origin[], IAtom atom) {
        AtomSite site = neighborSiteManager.getSite(atom);
        if (site == null) return;
        RectangularLattice lattice = neighborSiteManager.getLattice();
        lattice.latticeIndex(site.getLatticeArrayIndex(),latticeIndex);

        //color central site red
//        ((MySite)lattice.site(iterator.centralSite)).color = Color.RED;
        int nx = 5; //number of planes to draw before moving down to another line of planes
        int k = latticeIndex.length < 3 ? 0 : latticeIndex[2];
        int ox = k % nx;       //set origin for drawing each lattice plane
        int oy = (k - ox)/nx;
        ox = origin[0] + ox*spinWidth*(lattice.getSize()[0]) + 1;
        oy = origin[1] + oy*spinWidth*(lattice.getSize()[1]) + 1;
        //draw lattice plane
        g.setColor(atom.getPosition().getX(0) > 0 ? Color.green : Color.white);
        g.fillRect(ox+latticeIndex[0]*spinWidth,oy+latticeIndex[1]*spinWidth,spinWidth,spinWidth);
//        g.setColor(Color.black);
//        g.drawRect(ox+latticeIndex[0]*spinWidth,oy+latticeIndex[1]*spinWidth,spinWidth,spinWidth);
    }

    private int spinWidth;
    private final int[] latticeIndex;
    private final NeighborSiteManager neighborSiteManager;
}
