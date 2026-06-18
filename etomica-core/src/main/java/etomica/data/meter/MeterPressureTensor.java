/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.types.DataTensor;
import etomica.potential.PotentialCallbackVirialTensor;
import etomica.potential.compute.PotentialCallback;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.units.dimensions.Pressure;

public class MeterPressureTensor implements IDataSourcePotential {

    protected final DataTag tag;
    protected final DataTensor data;
    protected final DataInfo dataInfo;
    protected final PotentialCompute potentialMaster;
    protected final Box box;
    protected final PotentialCallbackVirialTensor pc;
    protected double temperature;
    protected boolean doNonEquilibrium;
    protected boolean callComputeAll = true;

    /**
     * Creates a meter that uses atomic velocities (like when doNonEquilibrium is set
     * to true) instead of including ideal gas contribution.
     */
    public MeterPressureTensor(PotentialCompute potentialMaster, Box box) {
        this(potentialMaster, box, Double.NaN);
        doNonEquilibrium = true;
    }

    public MeterPressureTensor(PotentialCompute potentialMaster, Box box, double temperature) {
        this.potentialMaster = potentialMaster;
        this.box = box;
        this.temperature = temperature;
        Space space = box.getSpace();
        data = new DataTensor(space);
        tag = new DataTag();
        dataInfo = new DataTensor.DataInfoTensor("Pressure", Pressure.dimension(space.D()), space);
        dataInfo.addTag(tag);
        pc = new PotentialCallbackVirialTensor(space);
    }

    public DataTag getTag() {
        return tag;
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public void setTemperature(double newT) {
        temperature = newT;
    }

    public void setDoNonEquilibrium(boolean doNonEquilibrium) {
        this.doNonEquilibrium = doNonEquilibrium;
    }

    /**
     * Computes total pressure in box by summing virial over all pairs, and adding
     * ideal-gas contribution.
     */
    public IData getData() {
        if (box == null) {
            throw new IllegalStateException("You must call setBox before using this class");
        }

        if (callComputeAll) {
            pc.reset();
            potentialMaster.computeAll(false, pc);
        }
        data.x.E(pc.getVirialTensor());

        IAtomList leafList = box.getLeafList();

        Tensor t = box.getSpace().makeTensor();
        if (doNonEquilibrium) {

            // use the velocities
            for (IAtom iAtom : leafList) {
                IAtomKinetic atom = (IAtomKinetic) iAtom;
                t.Ev1v2(atom.getVelocity(), atom.getVelocity());
                t.TE(atom.getType().getMass());
                data.x.PE(t);
            }
        } else {
            // or just include ideal gas term
            Vector v = box.getSpace().makeVector();
            v.E(1);
            t.diagE(v);
            data.x.PEa1Tt1(leafList.size() * temperature, t);
        }

        data.x.TE(1 / box.getBoundary().volume());
        return data;
    }

    public boolean needsPairCallback() {
        return true;
    }

    public PotentialCallback getPotentialCallback() {
        return pc;
    }

    @Override
    public void doCallComputeAll(boolean callComputeAll) {
        this.callComputeAll = callComputeAll;
    }
}
