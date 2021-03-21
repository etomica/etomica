/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.atom.IAtom;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.integrator.IntegratorBoxFasterer;
import etomica.molecule.IMolecule;
import etomica.potential.P2HardGeneric;
import etomica.units.dimensions.Null;

/**
 * Calculates the Boltzmann factor at a position within a box a molecule of a
 * particular species would have if it existed at that point.  Only works for
 * monatomic molecules (no rotation is attempted).
 *
 * @author Andrew Schultz
 */
public class DataSourcePositionedBoltzmannFactorFasterer implements DataSourceMolecular {

    protected final DataInfoDouble dataInfo;
    protected final DataDouble data;
    protected IntegratorBoxFasterer integrator;
    protected final DataTag tag;
    protected final P2HardGeneric p2Ghost, p2GhostHead, p2GhostTail;
    protected double headEpsilon;

    public DataSourcePositionedBoltzmannFactorFasterer(IntegratorBoxFasterer integrator, P2HardGeneric p2Ghost,
                                                       P2HardGeneric p2GhostHead, P2HardGeneric p2GhostTail, double headEpsilon) {
        data = new DataDouble();
        dataInfo = new DataInfoDouble("chemical potential", Null.DIMENSION);
        tag = new DataTag();
        this.integrator = integrator;
        this.p2Ghost = p2Ghost;
        this.p2GhostHead = p2GhostHead;
        this.p2GhostTail = p2GhostTail;
        setHeadEpsilon(headEpsilon);
    }

    public void setHeadEpsilon(double newHeadEpsilon) {
        headEpsilon = newHeadEpsilon;
    }

    public IData getData(IMolecule molecule) {
        p2Ghost.setEnergyForState(0, Double.POSITIVE_INFINITY);
        p2Ghost.setEnergyForState(1, -1.0);
        p2GhostHead.setEnergyForState(0, Double.POSITIVE_INFINITY);
        p2GhostHead.setEnergyForState(1, -1.0);
        p2GhostTail.setEnergyForState(0, Double.POSITIVE_INFINITY);

        IAtom atom = molecule.getChildList().get(0);
        double u = integrator.getPotentialCompute().computeOne(atom);
        double temp = integrator.getTemperature();
        data.x = Math.exp(-u / temp);

        p2Ghost.setEnergyForState(0, 0);
        p2Ghost.setEnergyForState(1, 0);
        p2GhostHead.setEnergyForState(0, 0);
        p2GhostHead.setEnergyForState(1, 0);
        p2GhostTail.setEnergyForState(0, 0);

        return data;
    }

    public IDataInfo getMoleculeDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

}
