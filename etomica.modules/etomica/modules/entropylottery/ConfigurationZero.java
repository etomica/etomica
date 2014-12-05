/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.entropylottery;

import etomica.action.MoleculeActionTranslateTo;
import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IVectorMutable;
import etomica.config.Configuration;
import etomica.space.ISpace;

/**
 * Configuration that simply puts all the Atoms at 0, or -0.5 if the box size
 * is even (which makes the entropy lottery binning happy).
 * @author Andrew Schultz
 */
public class ConfigurationZero implements Configuration, java.io.Serializable {

	private final ISpace space;

    public ConfigurationZero(ISpace _space) {
        super();
        this.space = _space;
    }

    public void initializeCoordinates(IBox box) {
        MoleculeActionTranslateTo atomActionTranslateTo = new MoleculeActionTranslateTo(space);
        IVectorMutable work = space.makeVector();
        work.E(0.0);
        int intD = (int)Math.round(box.getBoundary().getBoxSize().getX(0));
        if (intD % 2 == 0) {
            work.E(-0.5);
        }
        atomActionTranslateTo.setDestination(work);

        IAtomList leafList = box.getLeafList();
        for (int i=0; i<leafList.getAtomCount(); i++) {
            leafList.getAtom(i).getPosition().E(work);
        }
    }

    private static final long serialVersionUID = 2L;
}
