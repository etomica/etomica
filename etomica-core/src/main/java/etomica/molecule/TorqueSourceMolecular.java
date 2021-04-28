/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.molecule;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Boundary;
import etomica.space.Vector;

public class TorqueSourceMolecular {

    protected final PotentialCompute potentialCompute;
    protected final Boundary boundary;
    protected final Vector dr;
    protected final Vector torque;

    public TorqueSourceMolecular(Box box, PotentialCompute potentialCompute) {
        this.potentialCompute = potentialCompute;
        this.boundary = box.getBoundary();
        dr = box.getSpace().makeVector();
        torque = box.getSpace().makeVector();
    }

    public Vector getTorque(IMolecule molecule) {
        torque.E(0);
        Vector com = MoleculePositionCOMPBC.com(boundary, molecule);
        Vector[] torques = potentialCompute.getTorques();
        Vector[] forces = potentialCompute.getForces();
        for (IAtom a : molecule.getChildList()) {
            if (torques != null) {
                torque.PE(torques[a.getLeafIndex()]);
            }
            dr.Ev1Mv2(a.getPosition(), com);
            boundary.nearestImage(dr);
            dr.XE(forces[a.getLeafIndex()]);
            torque.ME(dr);
        }
        return torque;
    }
}
