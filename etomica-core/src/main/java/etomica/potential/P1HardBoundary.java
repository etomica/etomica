/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.space.Boundary;
import etomica.box.Box;
import etomica.space.Vector;
import etomica.graphics.Drawable;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Length;
import etomica.util.Debug;

/**
 * Potential that places hard repulsive walls coinciding with the
 * boundary of the box, which is assumed to be rectangular in shape.
 *
 * @author David Kofke
 */
public class P1HardBoundary implements PotentialHard, Drawable {
    
    private static final long serialVersionUID = 1L;
    private double collisionRadius = 0.0;
    private final Vector work;
    private int[] pixPosition;
    private int[] thickness;
    private int nominalThickness = 1;
    private boolean ignoreOverlap;
    private double lastVirial;
    private int lastCollisionDim;
    private final Tensor lastVirialTensor;
    protected Boundary boundary;
    
    public P1HardBoundary(Space space) {
        this(space, false);
    }
    
    public P1HardBoundary(Space space, boolean ignoreOverlap) {
        this.ignoreOverlap = ignoreOverlap;
        work = space.makeVector();
        lastVirialTensor = space.makeTensor();
        isActiveDim = new boolean[space.D()][2];
        for (int i=0; i<isActiveDim.length; i++) {
            isActiveDim[i][0] = true;
            isActiveDim[i][1] = true;
        }
    }

    public int nBody() {
        return 1;
    }

    public double getRange() {
        return collisionRadius;
    }

    public void setBox(Box box) {
        boundary = box.getBoundary();
    }

    public double energy(IAtomList a) {
        Vector dimensions = boundary.getBoxSize();
        Vector pos = a.getAtom(0).getPosition();
        for (int i=0; i<work.getD(); i++) {
            if (!isActiveDim[i][1]) {
                continue;
            }
            double rx = pos.getX(i);
            double dxHalf = 0.5*dimensions.getX(i);
            if((rx < -dxHalf+collisionRadius) || (rx > dxHalf-collisionRadius)) {
                 return Double.POSITIVE_INFINITY;
            }
        }
        return 0;
    }
     
    public double energyChange() {return 0.0;}
    
    public double collisionTime(IAtomList a, double falseTime) {
        IAtomKinetic atom = (IAtomKinetic)a.getAtom(0);
        work.E(atom.getPosition());
        Vector v = atom.getVelocity();
        work.PEa1Tv1(falseTime,v);
        Vector dimensions = boundary.getBoxSize();
        double tmin = Double.POSITIVE_INFINITY;
        for(int i=work.getD()-1; i>=0; i--) {
            double vx = v.getX(i);
            if(vx == 0.0) continue;
            double rx = work.getX(i);
            double dxHalf = 0.5*dimensions.getX(i);
            double t=0;
            if (vx > 0.0) {
                if (isActiveDim[i][1]) {
                    t = (dxHalf - rx - collisionRadius)/vx;
                }
                else {
                    continue;
                }
            }
            else if (isActiveDim[i][0]) {
                t = (-dxHalf -rx + collisionRadius)/vx;
            }
            else {
                continue;
            }
            if(t < tmin) tmin = t;
        }
        if (ignoreOverlap && tmin<0.0) tmin = 0.0;
        if (Debug.ON && tmin < 0.0) {
            System.out.println("t "+tmin+" "+atom+" "+work+" "+v+" "+boundary.getBoxSize());
            throw new RuntimeException("you screwed up");
        }
        return tmin + falseTime;
    }

    public void bump(IAtomList a, double falseTime) {
        IAtomKinetic atom = (IAtomKinetic)a.getAtom(0);
        work.E(atom.getPosition());
        Vector v = atom.getVelocity();
        work.PEa1Tv1(falseTime,v);
        Vector dimensions = boundary.getBoxSize();
        double delmin = Double.MAX_VALUE;
        int imin = 0;
        //figure out which component is colliding
        for(int i=work.getD()-1; i>=0; i--) {
            double rx = work.getX(i);
            double vx = v.getX(i);
            double dxHalf = 0.5*dimensions.getX(i);
            double del = (vx > 0.0) ? Math.abs(dxHalf - rx - collisionRadius) : Math.abs(-dxHalf - rx + collisionRadius);
            if(del < delmin) {
                delmin = del;
                imin = i;
            }
        }
        if (Debug.ON && collisionRadius > 0 && Math.abs(work.getX(imin)-collisionRadius+dimensions.getX(imin)*0.5)/collisionRadius > 1.e-9 
                && Math.abs(0.5*dimensions.getX(imin)-work.getX(imin)-collisionRadius)/collisionRadius > 1.e-9) {
            System.out.println(atom+" "+work+" "+dimensions);
            System.out.println("stop that");
        }
        v.setX(imin,-v.getX(imin));
        // dv = 2*NewVelocity
        double newP = atom.getPosition().getX(imin) - falseTime*v.getX(imin)*2.0;
        atom.getPosition().setX(imin,newP);
        double dp = 2.0/(atom.getType().rm())*(-v.getX(imin));
        lastVirial = dp;
        lastCollisionDim = imin;
    }//end of bump
    
    public double lastCollisionVirial() {
        // return 0 because the wall is not a molecule!
        return 0;
    }
    
    public Tensor lastCollisionVirialTensor() {
        // let's hope people only call this on purpose.  It should really be 0.
        // we're really returning the change in momentum 
        lastVirialTensor.E(0);
        lastVirialTensor.setComponent(lastCollisionDim, lastCollisionDim, lastVirial);
        return lastVirialTensor;
    }

    public double lastWallVirial() {
        double area = 1.0;
        final Vector dimensions = boundary.getBoxSize();
        for (int i=0; i<dimensions.getD(); i++) {
            if (i != lastCollisionDim) {
                area *= (dimensions.getX(i)-collisionRadius*2.0);
            }
        }
        double s = lastVirial / area;
        return s;
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
    public double getCollisionRadius() {return collisionRadius;}
    /**
     * Indicates collision radius has dimensions of Length.
     */
    public Dimension getCollisionRadiusDimension() {return Length.DIMENSION;}

    public void setActive(int dim, boolean first, boolean isActive) {
        isActiveDim[dim][first?0:1] = isActive;
    }
    
    public void setLongWall(int dim, boolean first, boolean longWall) {
        if (longWallDim == null) {
            longWallDim = new boolean[2][2];
            pixPosition = new int[2];
            thickness = new int[2];
        }
        longWallDim[dim][first?0:1] = longWall;
    }
    
    public void setDrawingThickness(int newThickness) {
        nominalThickness = newThickness;
    }

    public void draw(java.awt.Graphics g, int[] origin, double toPixel) {
        if (boundary == null) return;
        g.setColor(java.awt.Color.gray);
        // if not 2D serious problems!
        for (int i=0; i<2; i++) {
            for (int j=0; j<2; j++) {
                if (!isActiveDim[i][j]) continue;
                thickness[i] = nominalThickness;
                pixPosition[i] = origin[i] - thickness[i]/2;
                if (longWallDim == null || !longWallDim[i][j]) {
                    pixPosition[1-i] = origin[1-i];
                    thickness[1-i] = (int)(boundary.getBoxSize().getX(1-i)*toPixel);
                }
                else {
                    pixPosition[1-i] = 0;
                    thickness[1-i] = Integer.MAX_VALUE;
                }
                if (j==1) {
                    pixPosition[i] += (int)(boundary.getBoxSize().getX(i)*toPixel)-1;
                }
                g.fillRect(pixPosition[0],pixPosition[1],thickness[0],thickness[1]);
            }
        }
    }
    
    private boolean[][] isActiveDim;
    private boolean[][] longWallDim;
}
