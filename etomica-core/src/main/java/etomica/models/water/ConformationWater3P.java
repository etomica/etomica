/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.water;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.config.IConformation;
import etomica.space.Space;

/**
 * Conformation for 3-point water molecule.
 */
public class ConformationWater3P implements IConformation {

    public static final double bondLengthOH = 1.0;
    public static final double angleHOH = 109.5*Math.PI/180.;

    public ConformationWater3P(Space space) {
        this.space = space;
    }
    
    public void initializePositions(IAtomList list){
        
        IAtom o = list.get(2);
        o.getPosition().E(new double[] {0, 0, 0.0});

        double x = bondLengthOH*Math.sin(0.5*angleHOH);
        double y = bondLengthOH*Math.cos(0.5*angleHOH);
        
        IAtom h1 = list.get(0);
        h1.getPosition().E(new double[] {-x, y, 0.0});
                
        IAtom h2 = list.get(1);
        h2.getPosition().E(new double[] {+x, y, 0.0});
    }
    
    private static final long serialVersionUID = 1L;
    protected final Space space;
}
