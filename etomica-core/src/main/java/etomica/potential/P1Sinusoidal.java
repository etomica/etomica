/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtom;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.CompoundDimension;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Energy;
import etomica.units.dimensions.Length;

/**
 * Potential in which attaches a harmonic spring between each affected atom and
 * the nearest boundary in each direction.
 * <p>
 * This class has not been used or checked for correctness.
 *
 * @author David Kofke
 */

public class P1Sinusoidal implements IPotential1 {

    private final Space space;
    private double a, k;
    private double pi = Math.PI;
    private Vector shift;
    private boolean isBCC;

    public P1Sinusoidal(Space space, double a, double k, boolean isBCC, Vector shift) {
        this.space = space;
        this.a = a;
        this.k = k;
        this.shift = shift;
        this.isBCC = isBCC;
    }

    public void setSpringConstant(double springConstant) {
        k = springConstant;
    }
    
    public double getSpringConstant() {
        return k;
    }

    public Dimension getX0Dimension() {
        return Length.DIMENSION;
    }

    public Dimension getSpringConstantDimension() {
        return new CompoundDimension(new Dimension[]{Energy.DIMENSION, Length.DIMENSION}, new double[]{1, -2});
    }

    public double u(IAtom atom) {
        double x = atom.getPosition().getX(0) - shift.getX(0);
        double y = atom.getPosition().getX(1) - shift.getX(1);
        double z = atom.getPosition().getX(2) - shift.getX(2);

        if (isBCC) {
            double sin1 = Math.sin(pi*(y+z)/a);
            double sin2 = Math.sin(pi*(x+z)/a);
            double sin3 = Math.sin(pi*(x+y)/a);
            double u = k*(sin1*sin1 + sin2*sin2 + sin3*sin3);
            return u;
        } else {
            double sin1 = Math.sin(pi*(-x+y+z)/a);
            double sin2 = Math.sin(pi*(x-y+z)/a);
            double sin3 = Math.sin(pi*(x+y-z)/a);
            double u = k*(sin1*sin1 + sin2*sin2 + sin3*sin3);
            return u;
        }
    }

    public double udu(IAtom atom, Vector f) {
        double x = atom.getPosition().getX(0) - shift.getX(0);
        double y = atom.getPosition().getX(1) - shift.getX(1);
        double z = atom.getPosition().getX(2) - shift.getX(2);

        double u;
        if (isBCC) {
            double sin1 = Math.sin(pi*(y+z)/a);
            double sin2 = Math.sin(pi*(x+z)/a);
            double sin3 = Math.sin(pi*(x+y)/a);
            u = k*(sin1*sin1 + sin2*sin2 + sin3*sin3);

            double sin11 = Math.sin(2*pi*(y+z)/a);
            double sin22 = Math.sin(2*pi*(x+z)/a);
            double sin33 = Math.sin(2*pi*(x+y)/a);
            double Fx = -k*pi/a*(sin22+sin33);
            double Fy = -k*pi/a*(sin11+sin33);
            double Fz = -k*pi/a*(sin11+sin22);
            Vector Fxyz = space.makeVector();
            Fxyz.setX(0, Fx);
            Fxyz.setX(1, Fy);
            Fxyz.setX(2, Fz);
            f.PE(Fxyz);

            return u;
        } else {
            double sin1 = Math.sin(pi*(-x+y+z)/a);
            double sin2 = Math.sin(pi*(x-y+z)/a);
            double sin3 = Math.sin(pi*(x+y-z)/a);
            u = k*(sin1*sin1 + sin2*sin2 + sin3*sin3);

            double sin11 = Math.sin(2*pi*(-x+y+z)/a);
            double sin22 = Math.sin(2*pi*(x-y+z)/a);
            double sin33 = Math.sin(2*pi*(x+y-z)/a);
            double Fx = -k*pi/a*(-sin11+sin22+sin33);
            double Fy = -k*pi/a*(sin11-sin22+sin33);
            double Fz = -k*pi/a*(sin11+sin22-sin33);
            Vector Fxyz = space.makeVector();
            Fxyz.setX(0, Fx);
            Fxyz.setX(1, Fy);
            Fxyz.setX(2, Fz);
            f.PE(Fxyz);

            return u;
        }
    }

}