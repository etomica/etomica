/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.config;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.space.Space;


public class ConformationWater implements IConformation, java.io.Serializable {

    private static final long serialVersionUID = 1L;
    private double bondLengthOH = 4.0;
    private double angleHOH = 109.5*Math.PI/180.;

    public ConformationWater(Space space) {
        this.space = space;
    }
    
    public void initializePositions(IAtomList list) {
        
        double x = 6.0;
        double y = 6.0;
        
        IAtom o = list.getAtom(0);
        o.getPosition().E(new double[] {x, y, 0.0});
               
        IAtom h1 = list.getAtom(1);
        h1.getPosition().E(new double[] {x+bondLengthOH, y, 0.0});
                
        IAtom h2 = list.getAtom(2);
        h2.getPosition().E(new double[] {x+bondLengthOH*Math.cos(angleHOH), y+bondLengthOH*Math.sin(angleHOH), 0.0});
    }
        
    protected final Space space;
}
