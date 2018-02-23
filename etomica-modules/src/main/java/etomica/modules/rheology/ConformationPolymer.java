/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.rheology;

import etomica.atom.IAtomList;
import etomica.util.random.IRandom;
import etomica.space.Vector;
import etomica.config.IConformation;
import etomica.space.Space;

/**
 * Places polymer in a blob near the origin.
 *
 * @author Andrew Schultz
 */
public class ConformationPolymer implements IConformation {
    
    public ConformationPolymer(Space space, IRandom random) {
        this.random = random;
        r = space.makeVector();
    }

    public void initializePositions(IAtomList atomList) {
        r.E(0);
        for (int i = 0; i<atomList.size(); i++) {
            Vector p = atomList.get(i).getPosition();
            for (int j=0; j<p.getD(); j++) {
                p.setX(j, random.nextGaussian());
            }
            r.PE(p);
        }
        r.TE(-1.0/atomList.size());
        for (int i = 0; i<atomList.size(); i++) {
            atomList.get(i).getPosition().PE(r);
        }
    }

    protected final IRandom random;
    protected final Vector r;
}
