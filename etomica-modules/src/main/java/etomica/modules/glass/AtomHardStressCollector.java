/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.glass;

import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.IntegratorHardFasterer;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;

public class AtomHardStressCollector implements IntegratorHardFasterer.CollisionListener, AtomStressSource {

    private boolean velIncluded = false;
    private final double[][][] stress;
    private final Tensor stressPair;
    private final TimeSource timeSource;
    private double t0;

    public AtomHardStressCollector(TimeSource timeSource) {
        this.timeSource = timeSource;
        Box box = timeSource.getBox();
        int n = box.getLeafList().size();
        int D = box.getSpace().D();
        stress = new double[n][D][D];
        stressPair = box.getSpace().makeTensor();
        t0 = timeSource.getCurrentTime();
    }

    public double[][][] getStress() {
        if (!velIncluded) {
            IAtomList atoms = timeSource.getBox().getLeafList();
            double t = timeSource.getCurrentTime();
            double idt = 1 / (t - t0);
            int D = stress[0].length;
            for (int i = 0; i < atoms.size(); i++) {
                IAtomKinetic a = (IAtomKinetic) atoms.get(i);
                Vector velocity = a.getVelocity();
                for (int j = 0; j < D; j++) {
                    for (int k = 0; k < D; k++) {
                        stress[i][j][k] *= -idt;
                    }
                }
                stressPair.Ev1v2(velocity, velocity);
                stressPair.TE(a.getType().getMass());
                for (int j = 0; j < D; j++) {
                    for (int k = 0; k < D; k++) {
                        stress[i][j][k] = stress[i][j][k] / 2 + stressPair.component(j, k);
                    }
                }
            }
            t0 = t;
            velIncluded = true;
        }

        return stress;
    }

    @Override
    public void pairCollision(IAtomKinetic atom1, IAtomKinetic atom2, Vector rij, Vector dv, double virial, double tCollision) {
        if (velIncluded) {
            final int l = stress[0].length;
            for (int i = 0; i < stress.length; i++) {
                for (int j = 0; j < l; j++) {
                    for (int k = 0; k < l; k++) {
                        stress[i][j][k] = 0;
                    }
                }
            }
            velIncluded = false;
        }
        Tensor t = Space.getInstance(atom1.getPosition().getD()).makeTensor();
        t.Ev1v2(rij, rij);
        t.TE(virial / rij.squared());
        int idx1 = atom1.getLeafIndex();
        int idx2 = atom2.getLeafIndex();
        int D = t.D();
        for (int i = 0; i < D; i++) {
            for (int j = 0; j < D; j++) {
                stress[idx1][i][j] += t.component(i, j);
                stress[idx2][i][j] += t.component(i, j);
            }
        }
    }
}
