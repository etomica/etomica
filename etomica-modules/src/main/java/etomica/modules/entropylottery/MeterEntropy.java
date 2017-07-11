/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.entropylottery;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.space.Vector;
import etomica.data.DataSourceScalar;
import etomica.data.IEtomicaDataSource;
import etomica.units.dimensions.Null;

/**
 * Calculates the entropy of the distribution of Atoms in a Box
 * using Stirling's approximation.
 * @author Andrew Schultz
 */
public class MeterEntropy extends DataSourceScalar implements IEtomicaDataSource {

    public MeterEntropy() {
        super("entropy", Null.DIMENSION);
        atomCount = new int[0];
    }

    public double getDataAsScalar() {
        Vector dimensions = box.getBoundary().getBoxSize();
        if (atomCount.length != (int)Math.round(dimensions.getX(0))) {
            atomCount = new int[(int)Math.round(dimensions.getX(0))];
        }
        for (int i=0; i<atomCount.length; i++) {
            atomCount[i] = 0;
        }
        IAtomList leafList = box.getLeafList();
        for (int i=0; i<leafList.getAtomCount(); i++) {
            IAtom a = leafList.getAtom(i);
            int x = (int)Math.round(a.getPosition().getX(0)+dimensions.getX(0)*0.5-0.5);
            atomCount[x]++;
        }
        
        double sum = 0;
        double atomTotal = box.getLeafList().getAtomCount();
        for (int i=0; i<atomCount.length; i++) {
            if (atomCount[i] > 0) {
                sum += atomCount[i] * Math.log(atomCount[i]/atomTotal);
            }
        }
        return -sum;
    }

    public Box getBox() {
        return box;
    }

    public void setBox(Box newBox) {
        box = newBox;

    }

    private static final long serialVersionUID = 1L;
    protected Box box;
    protected int[] atomCount;
}
