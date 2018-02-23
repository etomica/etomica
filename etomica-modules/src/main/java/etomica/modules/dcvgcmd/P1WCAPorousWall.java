/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.dcvgcmd;

import etomica.atom.IAtomList;
import etomica.space.Vector;
import etomica.potential.Potential1;
import etomica.potential.PotentialSoft;
import etomica.space.Space;
import etomica.space.Tensor;

/**
 * This acts as a 1-body WCA potential wall perpendicular to the z direction 
 * with a hole.  The size and placement of the hole can be set.  The interior
 * of the hole is discontinuous; an atom in the hole that attempts to move out 
 * the hole and into the wall will experience a discontinuity. 
 */
public class P1WCAPorousWall extends Potential1 implements PotentialSoft {

    private static final long serialVersionUID = 1L;
    private final Vector[] gradient;
    private double sigma, sigma2;
    private double epsilon;
    private double cutoff, cutoff2;
    private double poreRadius, poreRadius2;
    private Vector[] poreCenters;
    private double z;

    public P1WCAPorousWall(Space space) {
        this(space, 1.0, 1.0);
    }

    public P1WCAPorousWall(Space space, double sigma, double epsilon) {
        super(space);
        setSigma(sigma);
        setEpsilon(epsilon);
        gradient = new Vector[1];
        gradient[0] = space.makeVector();
    }

    public double getRange() {
        return cutoff;
    }

    public double energy(IAtomList atom) {
        Vector r = atom.get(0).getPosition();
        double rz = r.getX(2);
        double dz2 = (z - rz);
        dz2 *= dz2;
        if(dz2 > cutoff2 || inPore(r)) return 0.0;
        return energy(dz2);
    }//end of energy

    private double energy(double r2) {
        r2 = sigma2 / r2;
        double r6 = r2 * r2 * r2;
        return 4 * epsilon * r6 * (r6 - 1.0) + epsilon;
    }
    
    public double virial(IAtomList atoms) {
        return 0.0;
    }
    
    private boolean inPore(Vector r) {
        for(int i=0; i<poreCenters.length; i++) {
            double dx = r.getX(0) - poreCenters[i].getX(0);
            double dy = r.getX(1) - poreCenters[i].getX(1);
            double r2 = dx*dx + dy*dy;
            if(r2 < poreRadius2) return true;
        }
        return false;
    }

    private double gradient(double r2) {
        r2 = sigma2 / r2;
        double r6 = r2 * r2 * r2;
        return -48 * epsilon * r6 * (r6 - 0.5);
    }

    public Vector[] gradient(IAtomList atom) {
        Vector r = atom.get(0).getPosition();
        double rz = r.getX(2);
        double dz2 = (z - rz);
        dz2 *= dz2;
        double gradz = 0.0;
        if(dz2 < cutoff2 && !inPore(r)) {
            gradz = gradient(dz2);
            if(z > rz) gradz = -gradz;
        }
        gradient[0].setX(2, gradz);
        return gradient;
    }
    
    public Vector[] gradient(IAtomList atom, Tensor pressureTensor) {
        return gradient(atom);
    }

    /**
     * Returns the radius.
     * 
     * @return double
     */
    public double getSigma() {
        return sigma;
    }

    /**
     * Sets the radius.
     * 
     * @param radius
     *            The radius to set
     */
    public void setSigma(double radius) {
        this.sigma = radius;
        sigma2 = radius * radius;
        cutoff = radius * Math.pow(2, 1. / 6.);
        cutoff2 = cutoff*cutoff;
    }

    /**
     * @return Returns the epsilon.
     */
    public double getEpsilon() {
        return epsilon;
    }

    /**
     * @param epsilon
     *            The epsilon to set.
     */
    public void setEpsilon(double epsilon) {
        this.epsilon = epsilon;
    }
    
    
    /**
     * @return Returns the poreRadius.
     */
    public double getPoreRadius() {
        return poreRadius;
    }
    /**
     * @param poreRadius The poreRadius to set.
     */
    public void setPoreRadius(double poreRadius) {
        this.poreRadius = poreRadius;
        poreRadius2 = poreRadius*poreRadius;
    }
    
    
    /**
     * @return Returns the z.
     */
    public double getZ() {
        return z;
    }
    /**
     * @param z The z to set.
     */
    public void setZ(double z) {
        this.z = z;
    }
    
    
    /**
     * @return Returns the poreCenters.
     */
    public Vector[] getPoreCenters() {
        return poreCenters;
    }
    /**
     * @param poreCenters The poreCenters to set.
     */
    public void setPoreCenters(Vector[] poreCenters) {
        this.poreCenters = poreCenters;
    }
}
