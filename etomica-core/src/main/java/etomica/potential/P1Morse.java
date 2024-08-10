/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtom;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.units.dimensions.CompoundDimension;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Energy;
import etomica.units.dimensions.Length;


public class P1Morse implements IPotential1 {

    private final Space space;

    private double re;
    private double epsilon;
    private double a;
    private Tensor eye;

    public P1Morse(Space space, double epsilon, double re, double a) {
        this.space = space;
        this.epsilon = epsilon;
        this.re = re;
        this.a = a;
        eye = space.makeTensor();
        eye.setComponent(0, 0, 1.0);
        eye.setComponent(1, 1, 1.0);
        eye.setComponent(2, 2, 1.0);

    }

    public Dimension getSpringConstantDimension() {
        return new CompoundDimension(new Dimension[]{Energy.DIMENSION, Length.DIMENSION}, new double[]{1, -2});
    }

    public double u(IAtom atom) {
        double r = Math.sqrt(atom.getPosition().squared());
        double expTerm = Math.exp(a*(re-r));
        return epsilon*(1-expTerm)*(1-expTerm);
    }

    public double udu(IAtom atom, Vector f) {
        double r = Math.sqrt(atom.getPosition().squared());
        double expTerm = Math.exp(a*(re-r));
        f.PEa1Tv1(-2*epsilon*a*expTerm*(1-expTerm)/r, atom.getPosition());
        return epsilon*(1-expTerm)*(1-expTerm);
    }

    @Override
    public Tensor d2u(IAtom atom) {
        Tensor Hii = space.makeTensor();
        double r = Math.sqrt(atom.getPosition().squared());
        double r2 = r*r;
        double r4 = r2*r2;
        double expTerm = Math.exp(a*(re-r));
        double du = r*2*epsilon*a*expTerm*(1-expTerm);  //rdu: r (2*epsilon*a*expTerm*(1-expTerm))
        double d2u = r2*2*epsilon*a*a*expTerm*(2*expTerm-1); //r2d2u: r^2 (k2 + 1/2 k4 r^2)
        Hii.Ev1v2(atom.getPosition(), atom.getPosition());
        Hii.TE((d2u-du)/r4);
        Hii.PEa1Tt1(du/r2, eye);
        return Hii;
    }

}