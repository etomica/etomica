/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.freeenergy.npath;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataSource;
import etomica.data.IEtomicaDataInfo;
import etomica.data.types.DataDoubleArray;
import etomica.meam.P2EAM;
import etomica.potential.IPotentialAtomic;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationEnergySum;
import etomica.potential.PotentialMaster;
import etomica.units.dimensions.Energy;

/**
 * Created by andrew on 4/12/17.
 */
public class DataSourceEnergies implements IDataSource {

    protected final DataDoubleArray data;
    protected final DataDoubleArray.DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
    protected final PotentialCalculationDUDW pcDUDW;
    protected final PotentialMaster potentialMaster;
    protected final IteratorDirective id;
    protected PotentialCalculationEnergies pc;
    protected Box box;

    public DataSourceEnergies(PotentialMaster potentialMaster) {
        this.potentialMaster = potentialMaster;
        data = new DataDoubleArray(3);
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("energies!", Energy.DIMENSION, new int[]{3});
        tag = new DataTag();
        dataInfo.addTag(tag);
        pc = new PotentialCalculationEnergies();
        pcDUDW = new PotentialCalculationDUDW();
        id = new IteratorDirective();
    }
    
    public void setPotentialCalculation(PotentialCalculationEnergies pc) {
        this.pc = pc;
    }

    public void setBox(Box box) {
        this.box = box;
    }

    @Override
    public IData getData() {
        pc.zeroSum();
        potentialMaster.calculate(box, id, pc);
        double[] x = data.getData();
        x[0] = pc.getSum1();
        x[1] = pc.getSum2();
        pcDUDW.reset();
        potentialMaster.calculate(box, id, pcDUDW);
        x[2] = pcDUDW.getSum();
        return data;
    }

    @Override
    public DataTag getTag() {
        return tag;
    }

    @Override
    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }

    public static class PotentialCalculationEnergies extends PotentialCalculationEnergySum {

        protected double sum1, sum2;

        public double getSum1() {
            return sum1;
        }

        public double getSum2() {
            return sum2;
        }

        public double getSum() {
            return sum1 + sum2;
        }

        public void zeroSum() {
            sum1 = sum2 = 0;
        }

        @Override
        public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
            double u = potential.energy(atoms);
            if (potential.nBody() == 1) {
                sum1 += u;
            }
            else {
                sum2 += u;
            }
        }
    }
    
    public static class PotentialCalculationEnergiesEAM extends PotentialCalculationEnergies {
        protected final P2EAM p2;
        public PotentialCalculationEnergiesEAM(P2EAM p2) {
            this.p2 = p2;
        }

        public void zeroSum() {
            super.zeroSum();
            p2.reset();
        }
        public double getSum2() {
            return sum2 + p2.energy1();
        }

        public double getSum() {
            return sum1 + sum2 + p2.energy1();
        }
    }
}
