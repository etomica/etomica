/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.zeno;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.space.Vector;
import etomica.util.random.IRandom;

public class MCMoveIdealChain extends MCMoveBox {

    protected final Box box;
    protected final IRandom random;
    protected double pBack, pForward;

    public MCMoveIdealChain(Box box, IRandom random, double pBack, double pForward) {
        this.box = box;
        this.random = random;
        this.pBack = pBack;
        this.pForward = pForward;
    }

    @Override
    public double energyChange() {
        return 0;
    }

    @Override
    public boolean doTrial() {
        IMoleculeList list = box.getMoleculeList();
        if (list.size() == 0) return false;
        int im = random.nextInt(list.size());
        IMolecule molecule = list.get(im);
        IAtomList atoms = molecule.getChildList();
        Vector x = box.getSpace().makeVector();
        Vector rPrev = atoms.get(0).getPosition();
        rPrev.E(x);
        int aPrev = random.nextInt(x.getD());
        x.setX(aPrev, random.nextInt(2)*2 - 1);
        for (int i=1; i<atoms.size(); i++) {
            Vector r = atoms.get(i).getPosition();
            r.Ev1Mv2(rPrev, x);
            double rnd = random.nextDouble();
            if (rnd < pBack) {
                x.TE(-1);
            }
            else if (rnd > pBack + pForward) {
                int a = aPrev;
                while (a == aPrev) a = random.nextInt(x.getD());
                x.E(0);
                x.setX(a, random.nextInt(2)*2 - 1);
                aPrev = a;
            }
            rPrev = r;
        }
        return false;
    }

    @Override
    public double getChi(double temperature) {
        return 1;
    }

    @Override
    public void acceptNotify() {

    }

    @Override
    public void rejectNotify() {

    }
}
