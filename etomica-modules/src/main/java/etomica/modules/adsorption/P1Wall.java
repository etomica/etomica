/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.adsorption;

import etomica.action.AtomActionRandomizeVelocity;
import etomica.atom.IAtomKinetic;
import etomica.integrator.IntegratorBox;
import etomica.potential.P1HardFieldGeneric;
import etomica.space.Vector;
import etomica.util.random.IRandom;

public class P1Wall extends P1HardFieldGeneric {

    public P1Wall(double L, double sigma, double range, double epsilon) {
        super(1, new double[3], new double[2]);
        setL(L);
        setSigma(sigma);
        setRange(range);
        setEpsilon(epsilon);
    }

    private void setPositions() {
        collisionPositions[0] = -L / 2 + sigma / 2;
        collisionPositions[1] = -L / 2 + sigma / 2 + range;
        collisionPositions[2] = +L / 2;
    }

    public void setL(double newL) {
        L = newL;
        setPositions();
    }

    public void setSigma(double newSigma) {
        sigma = newSigma;
        setPositions();
    }

    public double getSigma() {
        return sigma;
    }

    public void setRange(double newRange) {
        range = newRange;
        setPositions();
    }

    public void setEpsilon(double newEpsilon) {
        energies[0] = -newEpsilon;
    }

    public double getEpsilon() {
        return -energies[0];
    }

    public int bump(IAtomKinetic atom, int oldState, Vector r, double falseTime, Vector deltaP, double[] du) {
        int newState = super.bump(atom, oldState, r, falseTime, deltaP, du);
        if (r.getX(fieldDimension) < (collisionPositions[0] + collisionPositions[1]) / 2 && random.nextDouble() < pThermalize) {
            Vector v = atom.getVelocity();
            Vector vOld = Vector.d(v.getD());
            vOld.E(v);
            randomizer.setTemperature(integrator.getTemperature());
            randomizer.actionPerformed(atom);
            v.setX(fieldDimension, Math.abs(v.getX(fieldDimension)));
            vOld.ME(v);
            atom.getPosition().PEa1Tv1(falseTime, vOld);
        }
        return newState;
    }

    public void setThermalize(IntegratorBox integrator, double pThermalize, IRandom random) {
        this.integrator = integrator;
        this.pThermalize = pThermalize;
        this.random = random;
        randomizer = new AtomActionRandomizeVelocity(integrator.getTemperature(), random);
    }

    protected double L, range, sigma;
    protected double pThermalize;
    protected IntegratorBox integrator;
    protected IRandom random;
    protected AtomActionRandomizeVelocity randomizer;
}
