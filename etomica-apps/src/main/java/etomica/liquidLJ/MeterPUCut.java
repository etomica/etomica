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
import etomica.units.dimensions.Null;

/**
 * Meter for evaluation of the soft-potential pressure in a box.
 * Requires that temperature be set in order to calculation ideal-gas
 * contribution to pressure; default is to use zero temperature, which
 * causes this contribution to be omitted.
 *
 * @author David Kofke
 */

public class MeterPUCut implements IDataSource {
    
    protected final DataDoubleArray data;
    protected final DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
    protected IteratorDirective iteratorDirective;
    protected final PotentialCalculationSumCutoff pc, pcDADv2;
    protected PotentialMaster potentialMaster, potentialMasterDADv2;
    protected double temperature;
    protected Box box;
    private final int dim;
    
    public MeterPUCut(Space space, double[] cutoffs) {
        data = new DataDoubleArray(new int[]{cutoffs.length,4});
        dataInfo = new DataInfoDoubleArray("PU", Null.DIMENSION, new int[]{cutoffs.length,4});
        tag = new DataTag();
        dataInfo.addTag(tag);
    	dim = space.D();
        iteratorDirective = new IteratorDirective();
        iteratorDirective.includeLrc = false;
        pc = new PotentialCalculationSumCutoff(space, cutoffs);
        pcDADv2 = new PotentialCalculationSumCutoff(space, cutoffs);
    }

    public void setPotentialMaster(PotentialMaster newPotentialMaster) {
        potentialMaster = newPotentialMaster;
    }
    
    public void setPotentialMasterDADv2(PotentialMaster newPotentialMasterDADv2) {
        this.potentialMasterDADv2 = newPotentialMasterDADv2;
    }

    public void setTemperature(double newTemperature) {
        temperature = newTemperature;
    }

    public void setBox(Box newBox) {
        box = newBox;
        pc.setBox(box);
        pcDADv2.setBox(box);
    }

    /**
     * Computes total pressure in box by summing virial over all pairs, and adding
     * ideal-gas contribution.
     */
    public IData getData() {
        if (potentialMaster == null || box == null) {
            throw new IllegalStateException("You must call setIntegrator before using this class");
        }
        pc.zeroSums();
        potentialMaster.calculate(box, iteratorDirective, pc);
        double[] uSum = pc.getUSums();
        double[] vSum = pc.getVSums();

        double[] uSumDADv2 = uSum;
        double[] vSumDADv2 = vSum;
        if (potentialMasterDADv2 != null) {
            pcDADv2.zeroSums();
            potentialMasterDADv2.calculate(box, iteratorDirective, pcDADv2);
            uSumDADv2 = pcDADv2.getUSums();
            vSumDADv2 = pcDADv2.getVSums();
        }

        double[] x = data.getData();
        int j = 0;
        for (int i=0; i<uSum.length; i++) {
            double vol = box.getBoundary().volume();
            int N = box.getMoleculeList().size();
            double density = N / vol;

            double P = density*temperature - vSum[i]/(vol*dim);
            double U = uSum[i]/N;
            x[j+0] = U;
            x[j+1] = P;
            x[j+2] = U/(4*Math.pow(density,4));
            // dbA/drho at constant Y
            // (Z - 4u/T)/density  --  for SS, Z-1 = 4u/T
            // dbA/dv2 at constant Y
            // dbA/dv2 = dbA/rho * (-rho^3/2)

            U = uSumDADv2[i]/N;
            double Pex = -vSumDADv2[i]/(vol*dim);
            x[j+3] = -(Pex/(temperature*density) - 4 * U / (temperature))*density*density/2;
            
            j+=4;
        }
        return data;
    }

    public DataTag getTag() {
        return tag;
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

}
