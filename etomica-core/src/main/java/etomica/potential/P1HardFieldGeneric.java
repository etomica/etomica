/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential;

import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.box.Box;
import etomica.space.Vector;

public class P1HardFieldGeneric implements IPotentialHardField {

    protected int fieldDimension;
    protected double[] collisionPositions;
    protected double[] energies;

    public P1HardFieldGeneric(int fieldDimension, double[] collisionPositions, double[] energies) {
        if (collisionPositions.length < 2) {
            throw new IllegalArgumentException("Need at least two positions");
        }
        if (collisionPositions.length != energies.length + 1) {
            throw new IllegalArgumentException("Need one more position than energy");
        }
        this.fieldDimension = fieldDimension;
        this.collisionPositions = collisionPositions;
        this.energies = energies;
    }

    public void setCollisionPosition(int i, double p) {
        collisionPositions[i] = p;
    }

    public double getCollisionPosition(int i) {
        return collisionPositions[i];
    }

    public void setEnergy(int i, double u) {
        energies[i] = u;
    }

    public double getEnergyForState(int i) {
        // energy assumed to be infinity outside of our known positions
        return i < 0 || i >= energies.length ? Double.POSITIVE_INFINITY : energies[i];
    }

    public int getState(IAtomKinetic atom) {
        double x = atom.getPosition().getX(fieldDimension);
        if (x < collisionPositions[0]) throw new RuntimeException("state unknown");
        for (int i = 1; i < collisionPositions.length; i++) {
            if (x < collisionPositions[i]) {
                return i - 1;
            }
        }
        throw new RuntimeException("state unknown");
    }

    public double u(IAtom atom) {
        return getEnergyForState(getState((IAtomKinetic) atom));
    }

    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    public void setBox(Box box) {
    }

    public int bump(IAtomKinetic atom, int oldState, Vector r, double falseTime, double[] du) {
        double v = atom.getVelocity().getX(fieldDimension);
        int newState = oldState + (v > 0 ? +1 : -1);
        double uJump = getEnergyForState(newState) - getEnergyForState(oldState);
        double m = atom.getType().getMass();
        double ke = 0.5 * v * v * m;
        double vNew;
        if (ke > uJump) {
            vNew = Math.signum(v) * Math.sqrt(2 * (ke - uJump) / m);
        } else {
            newState = oldState;
            vNew = -v;
        }
        atom.getVelocity().setX(fieldDimension, vNew);
        Vector p = atom.getPosition();
        double x = p.getX(fieldDimension);
        p.setX(fieldDimension, x + (v - vNew) * falseTime);
        return newState;
    }

    public double collisionTime(IAtomKinetic atom, Vector r, Vector v, int state) {
        double vx = v.getX(fieldDimension);
        int i = state + (vx > 0 ? +1 : 0);
        double wallx = collisionPositions[i];
        return (wallx - r.getX(fieldDimension)) / vx;
    }
}
