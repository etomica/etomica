/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtom;
import etomica.space.Space;
import etomica.space.Vector;

public class P1Sinusoidal implements IPotential1 {

    private final Space space;
    private double a, k;
    private double pi = Math.PI;
    private Vector shift;
    private String strc;

    public P1Sinusoidal(Space space, double a, double k, String strc, Vector shift) {
        this.space = space;
        this.a = a;
        this.k = k;
        this.shift = shift;
        this.strc = strc;
    }

    public double u(IAtom atom) {
        double x = atom.getPosition().getX(0) - shift.getX(0);
        double y = atom.getPosition().getX(1) - shift.getX(1);
        double z = atom.getPosition().getX(2) - shift.getX(2);

        double sin1, sin2, sin3, sin4, sin5, sin6;
        double u = Double.NaN;
        if (strc == "SC") {
            sin1 = Math.sin(pi*x/a);
            sin2 = Math.sin(pi*y/a);
            sin3 = Math.sin(pi*z/a);
            u = k*(sin1*sin1+sin2*sin2+sin3*sin3);
        } else if (strc == "BCC") {
            sin1 = Math.sin(pi*(y+z)/a);
            sin2 = Math.sin(pi*(-y+z)/a);
            sin3 = Math.sin(pi*(x+z)/a);
            sin4 = Math.sin(pi*(-x+z)/a);
            sin5 = Math.sin(pi*(x+y)/a);
            sin6 = Math.sin(pi*(-x+y)/a);
            u = k*(sin1*sin1+sin2*sin2+sin3*sin3+sin4*sin4+sin5*sin5+sin6*sin6);
        } else if (strc == "FCC") {
            sin1 = Math.sin(pi*(x+y+z)/a);
            sin2 = Math.sin(pi*(-x+y+z)/a);
            sin3 = Math.sin(pi*(x-y+z)/a);
            sin4 = Math.sin(pi*(x+y-z)/a);
            u = k*(sin1*sin1+sin2*sin2+sin3*sin3+sin4*sin4);
        }

        return u;
    }

    public double udu(IAtom atom, Vector f) {
        double x = atom.getPosition().getX(0) - shift.getX(0);
        double y = atom.getPosition().getX(1) - shift.getX(1);
        double z = atom.getPosition().getX(2) - shift.getX(2);

        Vector Fxyz = space.makeVector();
        double sin1, sin2, sin3, sin4, sin5, sin6;
        double sin11, sin22, sin33, sin44, sin55, sin66;
        double Fx, Fy, Fz;
        double u = Double.NaN;
        if (strc == "SC") {
            sin1 = Math.sin(pi * x / a);
            sin2 = Math.sin(pi * y / a);
            sin3 = Math.sin(pi * z / a);
            u = k * (sin1 * sin1 + sin2 * sin2 + sin3 * sin3);

            sin11 = Math.sin(2 * pi * x / a);
            sin22 = Math.sin(2 * pi * y / a);
            sin33 = Math.sin(2 * pi * z / a);
            Fx = -k * pi / a * sin11;
            Fy = -k * pi / a * sin22;
            Fz = -k * pi / a * sin33;
            Fxyz.setX(0, Fx);
            Fxyz.setX(1, Fy);
            Fxyz.setX(2, Fz);
            f.PE(Fxyz);
        } else if (strc == "BCC") {
            sin1 = Math.sin(pi * (y + z) / a);
            sin2 = Math.sin(pi * (-y + z) / a);
            sin3 = Math.sin(pi * (x + z) / a);
            sin4 = Math.sin(pi * (-x + z) / a);
            sin5 = Math.sin(pi * (x + y) / a);
            sin6 = Math.sin(pi * (-x + y) / a);
            u = k * (sin1 * sin1 + sin2 * sin2 + sin3 * sin3 + sin4 * sin4 + sin5 * sin5 + sin6 * sin6);

            sin11 = Math.sin(2 * pi * (y + z) / a);
            sin22 = Math.sin(2 * pi * (-y + z) / a);
            sin33 = Math.sin(2 * pi * (x + z) / a);
            sin44 = Math.sin(2 * pi * (-x + z) / a);
            sin55 = Math.sin(2 * pi * (x + y) / a);
            sin66 = Math.sin(2 * pi * (-x + y) / a);
            Fx = -k * pi / a * (sin33 - sin44 + sin55 - sin66);
            Fy = -k * pi / a * (sin11 - sin22 + sin55 + sin66);
            Fz = -k * pi / a * (sin11 + sin22 + sin33 - sin44);
            Fxyz.setX(0, Fx);
            Fxyz.setX(1, Fy);
            Fxyz.setX(2, Fz);
            f.PE(Fxyz);
        } else if (strc == "FCC") {
            sin1 = Math.sin(pi * (x + y + z) / a);
            sin2 = Math.sin(pi * (-x + y + z) / a);
            sin3 = Math.sin(pi * (x - y + z) / a);
            sin4 = Math.sin(pi * (x + y - z) / a);
            u = k * (sin1 * sin1 + sin2 * sin2 + sin3 * sin3 + sin4 * sin4);

            sin11 = Math.sin(2*pi*(x+y+z)/a);
            sin22 = Math.sin(2*pi*(-x+y+z)/a);
            sin33 = Math.sin(2*pi*(x-y+z)/a);
            sin44 = Math.sin(2*pi*(x+y-z)/a);
            Fx = -k*pi/a*(sin11-sin22+sin33+sin44);
            Fy = -k*pi/a*(sin11+sin22-sin33+sin44);
            Fz = -k*pi/a*(sin11+sin22+sin33-sin44);
            Fxyz.setX(0, Fx);
            Fxyz.setX(1, Fy);
            Fxyz.setX(2, Fz);
            f.PE(Fxyz);
        }

        return u;
    }
}
