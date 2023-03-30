/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;

import etomica.action.controller.Controller;
import etomica.atom.IAtom;
import etomica.atom.IAtomOriented;
import etomica.space.Vector;

import java.awt.*;


public class DisplayBoxSpin2D extends DisplayBoxCanvas2D {

    /**
     * @param _box
     */
    public DisplayBoxSpin2D(DisplayBox _box, Controller controller) {
        super(_box, controller);
        spinWidth = 5;
    }
    
    protected void drawAtom(Graphics g, int origin[], IAtom atom) {
        Vector L = displayBox.getBox().getBoundary().getBoxSize();
        //color central site red
        int ox = origin[0] + 1;
        int oy = origin[1] + 1;
        //draw lattice plane
        g.setColor(((IAtomOriented)atom).getOrientation().getDirection().getX(0) > 0 ? Color.green : Color.white);
        double x = atom.getPosition().getX(0);
        double y = atom.getPosition().getX(1);
        g.fillRect(ox+(int)(x*spinWidth),oy+(int)(y*spinWidth),spinWidth,spinWidth);
//        g.setColor(Color.black);
//        g.drawRect(ox+latticeIndex[0]*spinWidth,oy+latticeIndex[1]*spinWidth,spinWidth,spinWidth);
    }

    private int spinWidth;
}
