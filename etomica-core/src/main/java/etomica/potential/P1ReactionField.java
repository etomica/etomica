/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential;

import etomica.atom.DipoleSourceAtomic;
import etomica.atom.IAtom;
import etomica.space.Vector;

public class P1ReactionField implements IPotential1 {

    protected final DipoleSourceAtomic dipoleSource;
    protected double epsilon;
    protected double cutoff;
    protected double fac;

    public P1ReactionField(DipoleSourceAtomic dipoleSource, double epsilon, double cutoff) {
        this.dipoleSource = dipoleSource;
        setDielectric(epsilon);
        setRange(cutoff);
    }

    public void setDielectric(double epsilon) {
        this.epsilon = epsilon;
        fac = 2 * (epsilon - 1) / (2 * epsilon + 1) / (cutoff * cutoff * cutoff);
    }

    public void setRange(double cutoff) {
        this.cutoff = cutoff;
        fac = 2 * (epsilon - 1) / (2 * epsilon + 1) / (cutoff * cutoff * cutoff);
    }

    @Override
    public double u(IAtom atom) {
        Vector iDipole = dipoleSource.getDipole(atom);
        return -0.5 * fac * iDipole.squared();
    }

    @Override
    public double udu(IAtom atom, Vector f) {
        return u(atom);
    }

    @Override
    public double uduTorque(IAtom atom, Vector f, Vector t) {
        return u(atom);
    }
}
