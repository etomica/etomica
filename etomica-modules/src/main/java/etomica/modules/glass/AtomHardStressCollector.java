/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.glass;

import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.IntegratorHard;
import etomica.space.Tensor;
import etomica.space.Vector;

public class AtomHardStressCollector implements IntegratorHard.CollisionListener, AtomStressSource {

    private boolean velIncluded = false;
    private final double[][] stress;
    private final Tensor stressPair;
    private final IntegratorHard integrator;
    private double t0;

    public AtomHardStressCollector(IntegratorHard integrator) {
        this.integrator = integrator;
        Box box = integrator.getBox();
        int n = box.getLeafList().size();
        int m = box.getSpace().D() == 3 ? 3 : 1;
        stress = new double[n][m];
        stressPair = box.getSpace().makeTensor();
        t0 = integrator.getCurrentTime();
    }

    public double[][] getStress() {
        if (!velIncluded) {
            IAtomList atoms = integrator.getBox().getLeafList();
            double t = integrator.getCurrentTime();
            double idt = 1 / (t - t0);
            for (int i = 0; i < atoms.size(); i++) {
                IAtomKinetic a = (IAtomKinetic) atoms.get(i);
                Vector velocity = a.getVelocity();
                for (int j = 0; j < stress[i].length; j++) {
                    stress[i][j] *= -idt;
                }
                stressPair.Ev1v2(velocity, velocity);
                stressPair.TE(a.getType().getMass());
                stress[i][0] += stressPair.component(0, 1);
                if (stress[i].length == 3) {
                    stress[i][1] += stressPair.component(0, 2);
                    stress[i][2] += stressPair.component(1, 2);
                }
            }
            t0 = t;
            velIncluded = true;
        }

        return stress;
    }

    public void collisionAction(IntegratorHard.Agent agent) {
        if (velIncluded) {
            final int l = stress[0].length;
            for (int i = 0; i < stress.length; i++) {
                for (int j = 0; j < l; j++) {
                    stress[i][j] = 0;
                }
            }
            velIncluded = false;
        }
        Tensor lcvt = agent.collisionPotential.lastCollisionVirialTensor();
        int idx0 = agent.atom.getLeafIndex();
        int idx1 = agent.collisionPartner.getLeafIndex();
        stress[idx0][0] += lcvt.component(0, 1);
        stress[idx1][0] += lcvt.component(0, 1);
        if (lcvt.D() == 3) {
            stress[idx0][1] += lcvt.component(0, 2);
            stress[idx1][1] += lcvt.component(0, 2);
            stress[idx0][2] += lcvt.component(1, 2);
            stress[idx1][2] += lcvt.component(1, 2);
        }
    }
}
