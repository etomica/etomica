/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.thermalconductivity;

import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.box.Box;
import etomica.potential.IPotential1;
import etomica.space.Vector;
import etomica.util.random.IRandom;

/**
 * Pair of thermal (diffuse) hard walls at the x faces of the box. The wall at
 * x = -Lx/2 is the hot wall and the wall at x = +Lx/2 is the cold wall. When a disk
 * collides with a wall its velocity is replaced by one drawn from the wall-flux
 * Maxwellian at the wall temperature: the normal component has density proportional
 * to v exp(-m v^2/2kT) and points into the fluid, and the tangential component is
 * Gaussian with variance kT/m (kB = 1).
 * <p>
 * The kinetic energy change of each collision is accumulated per wall (positive
 * means energy delivered to the fluid).
 */
public class P1HardThermalWall implements IPotential1 {

    protected final Box box;
    protected final IRandom random;
    protected double tHot, tCold;
    protected double radius;

    /** Net kinetic energy delivered to the fluid by the hot wall since last reset. */
    public double heatHot;
    /** Net kinetic energy delivered to the fluid by the cold wall (negative: energy removed). */
    public double heatCold;
    public long nCollHot, nCollCold;

    public P1HardThermalWall(Box box, IRandom random, double tHot, double tCold, double radius) {
        this.box = box;
        this.random = random;
        this.tHot = tHot;
        this.tCold = tCold;
        this.radius = radius;
    }

    public void setTemperatures(double tHot, double tCold) {
        this.tHot = tHot;
        this.tCold = tCold;
    }

    public double getTHot() {
        return tHot;
    }

    public double getTCold() {
        return tCold;
    }

    public void resetHeat() {
        heatHot = heatCold = 0;
        nCollHot = nCollCold = 0;
    }

    @Override
    public double u(IAtom atom) {
        double dxHalf = 0.5 * box.getBoundary().getBoxSize().getX(0);
        double x = atom.getPosition().getX(0);
        return (x < -dxHalf + radius || x > dxHalf - radius) ? Double.POSITIVE_INFINITY : 0;
    }

    @Override
    public double udu(IAtom atom, Vector f) {
        throw new RuntimeException("Can't compute force with hard potential");
    }

    @Override
    public int getState(IAtomKinetic atom) {
        return 0;
    }

    @Override
    public double getEnergyForState(int state) {
        return 0;
    }

    @Override
    public double collisionTime(IAtomKinetic atom, Vector r, Vector v, int state, double falseTime) {
        double vx = v.getX(0);
        if (vx == 0) return Double.POSITIVE_INFINITY;
        double dxHalf = 0.5 * box.getBoundary().getBoxSize().getX(0);
        double t = vx > 0 ? (dxHalf - radius - r.getX(0)) / vx : (-dxHalf + radius - r.getX(0)) / vx;
        // guard against tiny negative values from roundoff
        return t < 0 ? 0 : t;
    }

    @Override
    public int bump(IAtomKinetic atom, int oldState, Vector r, double falseTime, Vector deltaP, double[] du) {
        Vector v = atom.getVelocity();
        double m = atom.getType().getMass();
        double vxOld = v.getX(0);
        double vyOld = v.getX(1);
        boolean hot = r.getX(0) < 0;
        double s = Math.sqrt((hot ? tHot : tCold) / m);
        // normal component: density ~ v exp(-v^2/2s^2)  (1 - U avoids log(0))
        double vn = s * Math.sqrt(-2 * Math.log(1 - random.nextDouble()));
        double vxNew = hot ? vn : -vn;
        double vyNew = s * random.nextGaussian();
        double dKE = 0.5 * m * (vxNew * vxNew + vyNew * vyNew - vxOld * vxOld - vyOld * vyOld);

        Vector pos = atom.getPosition();
        pos.setX(0, pos.getX(0) + (vxOld - vxNew) * falseTime);
        pos.setX(1, pos.getX(1) + (vyOld - vyNew) * falseTime);
        v.setX(0, vxNew);
        v.setX(1, vyNew);

        deltaP.E(0);
        deltaP.setX(0, m * (vxNew - vxOld));
        deltaP.setX(1, m * (vyNew - vyOld));

        if (hot) {
            heatHot += dKE;
            nCollHot++;
        } else {
            heatCold += dKE;
            nCollCold++;
        }
        // integrator subtracts du from its kinetic energy; the walls are an energy reservoir
        du[0] = -dKE;
        return 0;
    }
}
