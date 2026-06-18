/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.pistoncylinder;

import etomica.atom.AtomSetSinglet;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.graphics.Drawable;
import etomica.potential.P1HardBoundary;
import etomica.space.Vector;
import etomica.units.dimensions.*;
import etomica.util.Debug;

/**
 * Potential that places hard repulsive walls that move and
 * accelerate subject to an external force field (pressure).
 */

public class P1HardMovingBoundary extends P1HardBoundary implements Drawable {

    /**
     * Constructor for a hard moving (and accelerating) boundary.
     *
     * @param wallDimension dimension which the wall is perpendicular to
     */
    public P1HardMovingBoundary(Box box, int wallDimension, double mass) {
        super(box.getSpace(), true, box);
        D = box.getSpace().D();
        wallD = wallDimension;
        setWallPosition(0);
        setMass(mass);
        force = 0.0;
    }

    public void setWallPosition(double p) {
        wallPosition = p;
    }

    public double getWallPosition() {
        return wallPosition;
    }

    public Dimension getWallPositionDimension() {
        return Length.DIMENSION;
    }

    public double getWallVelocity() {
        return wallVelocity;
    }

    public void setWallVelocity(double v) {
        wallVelocity = v;
    }

    public Dimension getWallVelocityDimension() {
        return new DimensionRatio("Velocity", Length.DIMENSION, Time.DIMENSION);
    }

    public void setPressure(double p) {
        pressure = p;
    }

    public double getPressure() {
        return pressure;
    }

    public Dimension getPressureDimension() {
        return Pressure.DIMENSION;
    }

    public void setStationary(boolean b) {
        if (b) {
            wallMass = Double.POSITIVE_INFINITY;
            wallVelocity = 0.0;
        } else {
            wallMass = setWallMass;
        }
    }

    public boolean isStationary() {
        return Double.isInfinite(wallMass);
    }

    /**
     * @return Returns the mass.
     */
    public double getMass() {
        return wallMass;
    }

    /**
     * @param mass The mass to set.
     */
    public void setMass(double mass) {
        wallMass = mass;
        setWallMass = mass;
    }

    public Dimension getMassDimension() {
        return Mass.DIMENSION;
    }

    public double energyChange() {
        return 0.0;
    }

    protected double getAcceleration() {
        double area = 1.0;
        final Vector dimensions = boundary.getBoxSize();
        for (int i = 0; i < D; i++) {
            if (i != wallD) {
                area *= (dimensions.getX(i) - collisionRadius * 2.0);
            }
        }
        force = pressure * area;

        return force / wallMass;   // atom acceleration - wall acceleration
    }

    @Override
    public double collisionTime(IAtomKinetic atom, Vector r, Vector v, int state, double falseTime) {
        double tboundary = super.collisionTime(atom, r, v, state, falseTime);

        double a = getAcceleration();
        double wallV = wallVelocity + a * falseTime;
        double wallP = wallPosition + wallVelocity * falseTime + 0.5 * a * falseTime * falseTime;

        double dr = r.getX(wallD) - wallP;
        double dv = atom.getVelocity().getX(wallD) - wallV;
        double da = 0 - a;

        if (Debug.ON && Debug.DEBUG_NOW && Debug.anyAtom(new AtomSetSinglet(atom))) {
            System.out.println(dr + " " + dv + " " + falseTime + " " + atom);
            System.out.println(atom.getVelocity().getX(wallD));
            System.out.println(atom.getPosition().getX(wallD));
        }
        double t = Double.POSITIVE_INFINITY;
        double discr = -1.0;
        if (dr * dv < 0.0 || dr * da < 0.0) {
            // either moving toward or accelerating toward each other
            if ((Debug.ON || ignoreOverlap) && Math.abs(dr) < collisionRadius && dr * dv < 0.0) {
                if (ignoreOverlap) return 0.001 * Math.abs(dr / dv);
                throw new RuntimeException("overlap " + atom + " " + dr + " " + dv + " " + da);
            }
            double drc;
            if (dr > 0.0) {
                drc = dr - collisionRadius;
            } else {
                drc = dr + collisionRadius;
            }
            discr = dv * dv - 2.0 * da * drc;
            if (discr >= 0.0) {
                discr = Math.sqrt(discr);
                if (dr * da < 0.0) {
                    t = -dv / da + discr / Math.abs(da);
                } else if (da == 0.0) {
                    if (dr * dv < 0.0) t = -drc / dv;
                } else if (dr * dv < 0.0 && dr * da > 0.0) {
                    t = -dv / da - discr / Math.abs(da);
                } else {
                    throw new RuntimeException("oops");
                }
            }
        }
        if (ignoreOverlap && t < 0.0) t = 0.001 * Math.abs(dr / dv);
        if (Debug.ON && (t < 0.0 || Debug.DEBUG_NOW && Debug.anyAtom(new AtomSetSinglet(atom)))) {
            System.out.println(atom + " " + da + " " + dr + " " + dv + " " + discr + " " + t + " " + (t + falseTime) + " " + (atom.getPosition().getX(wallD) + atom.getVelocity().getX(wallD) * (t + falseTime)) + " " + (wallPosition + wallVelocity * (t + falseTime) - 0.5 * da * (t + falseTime) * (t + falseTime)));
            if (t < 0) throw new RuntimeException("foo");
        }
        return Math.min(t, tboundary);
    }

    public double collisionTime(IAtomList atoms, double falseTime) {
        throw new RuntimeException("nope");
    }

    @Override
    public int bump(IAtomKinetic atom, int oldState, Vector r, double falseTime, Vector deltaP, double[] du) {
        double a = getAcceleration();
        double wallP = wallPosition + wallVelocity * falseTime + 0.5 * a * falseTime * falseTime;
        double x = r.getX(wallD);
        double dx = Math.abs(x - wallP) - collisionRadius;
        if (Math.abs(dx) > 1e-10) {
            return super.bump(atom, oldState, r, falseTime, deltaP, du);
        }

        double trueWallVelocity = wallVelocity + a * falseTime;
        Vector v = atom.getVelocity();
        double dp = 2.0 / (1 / wallMass + atom.getType().rm()) * (trueWallVelocity - v.getX(wallD));
        v.setX(wallD, v.getX(wallD) + dp * atom.getType().rm());
        atom.getPosition().setX(wallD, x - v.getX(wallD) * falseTime);
        wallVelocity -= dp / wallMass;
        wallPosition += dp / wallMass * falseTime;
        return 0;
    }

    public void bump(IAtomList atoms, double falseTime) {
        throw new RuntimeException("nope");
    }

    /**
     * Distance from the center of the sphere to the boundary at collision.
     */
    public void setCollisionRadius(double d) {
        if (d < 0) {
            throw new IllegalArgumentException("collision radius must not be negative");
        }
        collisionRadius = d;
    }

    /**
     * Distance from the center of the sphere to the boundary at collision.
     */
    public double getCollisionRadius() {
        return collisionRadius;
    }

    /**
     * Indicates collision radius has dimensions of Length.
     */
    public Dimension getCollisionRadiusDimension() {
        return Length.DIMENSION;
    }

    public void advanceAcrossTimeStep(double tStep) {
        double a = getAcceleration();
        wallPosition += wallVelocity * tStep + 0.5 * tStep * tStep * a;
        wallVelocity += tStep * a;
//        System.out.println("pressure => velocity "+a+" "+wallVelocity+" "+wallPosition+" "+tStep);
    }

    public void setThickness(double t) {
        thickness = t;
    }

    public void draw(java.awt.Graphics g, int[] origin, double toPixel) {
        super.draw(g, origin, toPixel);
        g.setColor(java.awt.Color.gray);
        double dx = boundary.getBoxSize().getX(0);
        double dy = boundary.getBoxSize().getX(1);
        int xP = origin[0] + (wallD == 0 ? (int) ((wallPosition + 0.5 * dx - thickness) * toPixel) : 0);
        int yP = origin[1] + (wallD == 1 ? (int) ((wallPosition + 0.5 * dy - thickness) * toPixel) : 0);
        int t = Math.max(1, (int) (thickness * toPixel));
        int wP = wallD == 0 ? t : (int) (toPixel * dx);
        int hP = wallD == 1 ? t : (int) (toPixel * dy);
        g.fillRect(xP, yP, wP, hP);
    }

    private final int D;
    private final int wallD;
    private double wallPosition;
    private double wallVelocity;
    private double wallMass;
    private double setWallMass;
    private double force;
    private double pressure;
    private double thickness = 0.0;
}
   
