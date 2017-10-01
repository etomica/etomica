/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.box.Box;
import etomica.data.*;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationHyperVirialSum;
import etomica.potential.PotentialCalculationVirialSum;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.units.dimensions.Null;

/**
 * Meter for measurement of the hypervirial, r^2 d2u/dr2 + r du/dr
 */
public class MeterHyperVirial implements IDataSource {

    protected final IDataInfo dataInfo;
    protected final IData data;
    protected final DataTag tag;
    protected final IteratorDirective iteratorDirective;
    protected final PotentialCalculationHyperVirialSum hyperVirial;
    protected PotentialMaster potentialMaster;
    protected Box box;

    private final PotentialCalculationVirialSum virial;
    protected double temperature;
    private final int dim;

    public MeterHyperVirial() {
        this(null, false);
    }

    /**
     * If full is true, then the pressure is also calculated and three values
     * are returned.  The first is the hypervirial value that would be returned
     * normally.  The second is the pressure.  The last is
     *   hypervirial/(dim*dim*n) + pressure/density
     * This is part of the expression for dP/drho 
     *
     * @param space
     * @param doFull
     */
    public MeterHyperVirial(Space space, boolean doFull) {
        iteratorDirective = new IteratorDirective();
        iteratorDirective.includeLrc = true;
        hyperVirial = new PotentialCalculationHyperVirialSum();
        tag = new DataTag();

        if (doFull) {
            dataInfo = new DataInfoDoubleArray("hyperVirial", Null.DIMENSION, new int[]{3});
            data = new DataDoubleArray(3);
            dim = space.D();
            virial = new PotentialCalculationVirialSum();
        }
        else {
            dataInfo = new DataInfoDouble("hyperVirial", Null.DIMENSION);
            data = new DataDouble();
            dim = 0;
            virial = null;
        }
        dataInfo.addTag(tag);
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    public void setPotentialMaster(PotentialMaster newPotentialMaster) {
        potentialMaster = newPotentialMaster;
    }
    
    public void setBox(Box newBox) {
        box = newBox;
    }

    public void setTemperature(double temperature) {
        this.temperature = temperature;
    }

    /**
     * Sets flag indicating whether calculated energy should include
     * long-range correction for potential truncation (true) or not (false).
     */
    public void setIncludeLrc(boolean b) {
    	iteratorDirective.includeLrc = b;
    }
    
    /**
     * Indicates whether calculated energy should include
     * long-range correction for potential truncation (true) or not (false).
     */
    public boolean isIncludeLrc() {
    	return iteratorDirective.includeLrc;
    }

	 /**
	  * Computes total pressure in box by summing virial over all pairs, and adding
	  * ideal-gas contribution.
	  */
    public IData getData() {
        hyperVirial.zeroSum();
        potentialMaster.calculate(box, iteratorDirective, hyperVirial);
        double x = hyperVirial.getSum();
        
        if (virial != null) {
            double[] y = ((DataDoubleArray)data).getData();
            y[0] = x;
            virial.zeroSum();
    
            potentialMaster.calculate(box, iteratorDirective, virial);
            //System.out.println("fac="+(1/(box.getBoundary().volume()*box.getSpace().D())));
            double V = box.getBoundary().volume();
            int n = box.getMoleculeList().getMoleculeCount();
            double rho = n/V; 
            double p = rho*temperature - virial.getSum()/(box.getBoundary().volume()*dim);
            y[1] = p;
            y[2] = x/(9*n) + p/rho;
        }
        else {
            ((DataDouble)data).x = x;
        }
        
        return data;
    }
}
