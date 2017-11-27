/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.liquidLJ;

import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataSource;
import etomica.data.IDataInfo;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.units.dimensions.Energy;

/**
 * Meter for evaluation of the potential energy in a box.
 * Includes several related methods for computing the potential energy of a single
 * atom or molecule with all neighboring atoms
 *
 * @author David Kofke
 */

public class MeterPotentialEnergyCutoff implements IDataSource {
    
    public MeterPotentialEnergyCutoff(PotentialMaster potentialMaster, Space space, double[] cutoffs) {
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("energy", Energy.DIMENSION, new int[]{cutoffs.length});
        tag = new DataTag();
        data = new DataDoubleArray(cutoffs.length);
        dataInfo.addTag(tag);
        iteratorDirective.includeLrc = false;
        potential = potentialMaster;
        iteratorDirective.setDirection(IteratorDirective.Direction.UP);
        energy = new PotentialCalculationEnergySumCutoff(space, cutoffs);
    }
    
    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

   /**
    * Computes total potential energy for box.
    * Currently, does not include long-range correction to truncation of energy
    */
    public IData getData() {
        if (box == null) throw new IllegalStateException("must call setBox before using meter");
        energy.setBox(box);
    	energy.zeroSum();
        potential.calculate(box, iteratorDirective, energy);
        System.arraycopy(energy.getSums(), 0, data.getData(), 0, dataInfo.getLength());
        return data;
    }
    /**
     * @return Returns the box.
     */
    public Box getBox() {
        return box;
    }
    /**
     * @param box The box to set.
     */
    public void setBox(Box box) {
        this.box = box;
    }

    protected Box box;
    protected final DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
    protected final DataDoubleArray data;
    protected final IteratorDirective iteratorDirective = new IteratorDirective();
    protected final PotentialCalculationEnergySumCutoff energy;
    protected final PotentialMaster potential;
}
