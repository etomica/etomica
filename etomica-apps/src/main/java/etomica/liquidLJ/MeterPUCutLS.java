/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.liquidLJ;
import etomica.box.Box;
import etomica.potential.PotentialMaster;
import etomica.atom.iterator.IteratorDirective;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
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
 
public class MeterPUCutLS implements IEtomicaDataSource {
    
    protected final DataDoubleArray data;
    protected final DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
    protected IteratorDirective iteratorDirective;
    protected final PotentialCalculationSumCutoffLS pc, pcDADv2;
    protected PotentialMaster potentialMaster, potentialMasterDADv2;
    protected double temperature;
    protected Box box;
    private final int dim;

    public MeterPUCutLS(Space space, int nCut) {
        data = new DataDoubleArray(new int[]{nCut,4});
        dataInfo = new DataInfoDoubleArray("PU", Null.DIMENSION, new int[]{nCut,4});
        tag = new DataTag();
        dataInfo.addTag(tag);
    	dim = space.D();
        iteratorDirective = new IteratorDirective();
        iteratorDirective.includeLrc = false;
        pc = new PotentialCalculationSumCutoffLS();
        pcDADv2 = new PotentialCalculationSumCutoffLS();
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
    }

    public IData getData() {
        if (potentialMaster == null || box == null) {
            throw new IllegalStateException("You must call setIntegrator before using this class");
        }
        pc.zeroSums();
        potentialMaster.calculate(box, iteratorDirective, pc);
        double[][] uvSums = pc.getSums();
        double[][] uvSumsDADv2 = uvSums;

        if (potentialMasterDADv2 != null) {
            pcDADv2.zeroSums();
            potentialMasterDADv2.calculate(box, iteratorDirective, pcDADv2);
            uvSumsDADv2 = pcDADv2.getSums();
        }

        double[] x = data.getData();
        int j = 0;
        for (int i=0; i<uvSums[0].length; i++) {
            double vol = box.getBoundary().volume();
            int N = box.getMoleculeList().getMoleculeCount();
            double density = N / vol;
            
            double P = density*temperature - uvSums[1][i]/(vol*dim);
            double U = uvSums[0][i]/N;
            x[j+0] = U;
            x[j+1] = P;
            x[j+2] = U/(4*Math.pow(density,4));

            // dbA/drho at constant Y
            // (Z - 4u/T)/density  --  for SS, Z-1 = 4u/T
            // dbA/dv2 at constant Y
            // dbA/dv2 = dbA/rho * (-rho^3/2)
            U = uvSumsDADv2[0][i]/N;
            P = -uvSumsDADv2[1][i]/(vol*dim);
            x[j+3] = -(P/(temperature*density) - 4 * U / (temperature))*density*density/2;
            
            j+=4;
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
