/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialCalculationPressureTensor;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.units.dimensions.Null;


/**
 * Meter that returns virial tensor components and also components of F dr
 */
public class MeterPressureStuff implements IDataSource {
    protected final DataTag tag;
    protected final DataInfoDoubleArray dataInfo;
    protected final DataDoubleArray data;
    protected final PotentialMaster potentialMaster;
    protected final IteratorDirective iteratorDirective;
    protected final PotentialCalculationForceSum pcForce;
    protected final PotentialCalculationPressureTensor pcPressure;
    protected final CoordinateDefinition coordinateDefinition;
    protected final AtomLeafAgentManager<Vector> forceManager;

    public MeterPressureStuff(Space space, PotentialMaster potentialMaster, CoordinateDefinition coordinateDefinition) {
        this.coordinateDefinition = coordinateDefinition;
        this.potentialMaster = potentialMaster;
        iteratorDirective = new IteratorDirective();
        iteratorDirective.includeLrc = false;
        tag = new DataTag();

        pcForce = new PotentialCalculationForceSum();
        forceManager = new AtomLeafAgentManager<>(a -> space.makeVector(), coordinateDefinition.getBox());
        pcForce.setAgentManager(forceManager);

        pcPressure = new PotentialCalculationPressureTensor(space);
        pcPressure.setTemperature(0);

        int n = 6;
        dataInfo = new DataInfoDoubleArray("Stuff", Null.DIMENSION, new int[]{n});
        dataInfo.addTag(tag);
        data = new DataDoubleArray(n);
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    public IData getData() {
        Box box = coordinateDefinition.getBox();

        double[] x = data.getData();

        pcPressure.zeroSum();
        potentialMaster.calculate(box, iteratorDirective, pcPressure);
        Tensor pTensor = pcPressure.getPressureTensor();

        // this thinks T=0, so it gives us the virial contribution (Wxx...)
        for (int i = 0; i < 3; i++) {
            x[i] = pTensor.component(i, i);
            x[3 + i] = 0;
        }

        // now F dr (Fx dx, Fy dy, Fz dz)
        pcForce.reset();
        potentialMaster.calculate(box, iteratorDirective, pcForce);
        IAtomList atoms = box.getLeafList();
        Vector dr0 = box.getSpace().makeVector();
        Vector r0 = atoms.get(0).getPosition();
        dr0.Ev1Mv2(r0, coordinateDefinition.getLatticePosition(atoms.get(0)));
        for (IAtom atom : atoms) {
            Vector f = forceManager.getAgent(atom);
            Vector r = atom.getPosition();
            Vector dr = box.getSpace().makeVector();
            dr.Ev1Mv2(r, coordinateDefinition.getLatticePosition(atom));
            dr.ME(dr0);
            box.getBoundary().nearestImage(dr);
            for (int i = 0; i < 3; i++) {
                x[3 + i] += f.getX(i) * dr.getX(i);
            }
        }

        // conventional
        // P = (Wxx + Wyy + Wzz) / 3V + T rho
        // real space
        // P = (Wxx + Wyy + Wzz) / 3V + T rho + Pqh + fV F dr

        return data;
    }

}