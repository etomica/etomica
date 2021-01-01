/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.sam;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListenerMD;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Vector;

public class SamIntegratorListener implements IntegratorListenerMD {

    public AtomType atomTypeFixedY = null;
    protected final PotentialCompute compute;
    protected final Box box;

    public SamIntegratorListener(PotentialCompute compute, Box box) {
        this.compute = compute;
        this.box = box;
    }

    @Override
    public void integratorInitialized(IntegratorEvent e) {
    }

    @Override
    public void integratorStepStarted(IntegratorEvent e) {
        if (atomTypeFixedY == null) return;
        // need to do this before each step in case thermostat kicked them
        IAtomList atoms = box.getLeafList();
        for (int i = 0; i < atoms.size(); i++) {
            IAtomKinetic a = (IAtomKinetic) atoms.get(i);
            if (a.getType() == atomTypeFixedY) a.getVelocity().setX(1, 0);
        }
    }

    @Override
    public void integratorStepFinished(IntegratorEvent e) {
    }

    @Override
    public void integratorForcePrecomputed(IntegratorEvent e) {
    }

    @Override
    public void integratorForceComputed(IntegratorEvent e) {
        Vector[] forces = compute.getForces();
        IAtomList atoms = box.getLeafList();
        for (int i = 0; i < atoms.size(); i++) {
            Vector f = forces[i];
            IAtom a = atoms.get(i);
            for (int j = 0; j < f.getD(); j++) {
                if (a.getType() == atomTypeFixedY) f.setX(1, 0);
                if (Math.abs(f.getX(i)) > 2.5e9) {
                    double fnew = 2.5e9 * Math.signum(f.getX(i));
                    f.setX(i, fnew);
                }
            }
        }
    }
}
