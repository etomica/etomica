/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.liquidLJ;

import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationEnergySum;
import etomica.potential.PotentialCalculationVirialSum;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.units.Null;

/**
 * Meter for evaluation of the soft-potential pressure in a box.
 * Requires that temperature be set in order to calculation ideal-gas
 * contribution to pressure; default is to use zero temperature, which
 * causes this contribution to be omitted.
 *
 * @author David Kofke
 */
 
public class MeterPU implements IEtomicaDataSource {
    
    protected final DataDoubleArray data;
    protected final DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
    protected IteratorDirective iteratorDirective;
    protected final PotentialCalculationVirialSum virial;
    protected final PotentialCalculationEnergySum energy;
    protected PotentialMaster potentialMaster;
    protected double temperature;
    protected Box box;
    private final int dim;
    
    public MeterPU(Space space) {
        data = new DataDoubleArray(4);
        dataInfo = new DataInfoDoubleArray("PU", Null.DIMENSION, new int[]{4});
        tag = new DataTag();
        dataInfo.addTag(tag);
    	dim = space.D();
        iteratorDirective = new IteratorDirective();
        iteratorDirective.includeLrc = true;
        virial = new PotentialCalculationVirialSum();
        energy = new PotentialCalculationEnergySum();
    }

    public void setPotentialMaster(PotentialMaster newPotentialMaster) {
        potentialMaster = newPotentialMaster;
    }
    
    public void setTemperature(double newTemperature) {
        temperature = newTemperature;
    }

    public void setBox(Box newBox) {
        box = newBox;
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
        if (potentialMaster == null || box == null) {
            throw new IllegalStateException("You must call setIntegrator before using this class");
        }
        virial.zeroSum();
        potentialMaster.calculate(box, iteratorDirective, virial);
        //System.out.println("fac="+(1/(box.getBoundary().volume()*box.getSpace().D())));
        double vol = box.getBoundary().volume();
        int N = box.getMoleculeList().getMoleculeCount();
        double density = N / vol;
        
        double P = density*temperature - virial.getSum()/(vol*dim);
        energy.zeroSum();
        potentialMaster.calculate(box, iteratorDirective, energy);
        double U = energy.getSum()/N;
        double[] x = data.getData();
        x[0] = U;
        x[1] = P;
        x[2] = U/(4*Math.pow(density,4));
        // dbA/drho at constant Y
        // (Z - 4u/T)/density  --  for SS, Z-1 = 4u/T
        // dbA/dv2 at constant Y
        // dbA/dv2 = dbA/rho * (-rho^3/2)
        x[3] = -(P/(temperature*density) - 1 - 4 * U / (temperature))*density*density/2;
        if (false) {
            // DADRY - DADY*Y*v^3
            // (Z - 4u/T)/rho  +  (4*u rho^4 /T^2) * Y * v^3
            // (Z - 4u/T)/rho  +  (4*u rho^4 /T^2) * (0.25 T / rho^4) * v^3
            // (Z - 4u/T)/rho  +  u/T * v^3
            // (Z - 4u/T + u v2 / T)/rho
            x[3] = (P/(temperature*density) - 1 - 4 * U / (temperature) - U/(density*density*temperature))/density;
        }
        return data;
    }

    public DataTag getTag() {
        return tag;
    }

    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }

}
