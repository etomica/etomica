/*
 * Created on Apr 9, 2004
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package etomica.modules.dcvgcmd;

import etomica.Default;
import etomica.EtomicaInfo;
import etomica.Space;
import etomica.atom.Atom;
import etomica.atom.AtomSet;
import etomica.potential.Potential1;
import etomica.potential.PotentialSoft;
import etomica.space.Vector;

/**
 * @author Owner
 * 
 * To change the template for this generated type comment go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */

public class P1WCAWall extends Potential1 implements PotentialSoft {

    private final Vector gradient;
    private double sigma;
    private double epsilon;
    private double cutoff;

    public P1WCAWall(Space space) {
        this(space, Default.ATOM_SIZE, Default.POTENTIAL_WELL);
    }

    public P1WCAWall(Space space, double sigma, double epsilon) {
        super(space);
        setSigma(sigma);
        setEpsilon(epsilon);
        gradient = space.makeVector();
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo(
                "WCA LJ Potential in the Z-Coordinate");
        return info;
    }

    public double energy(AtomSet atom) {
        Atom a = (Atom) atom;
        Vector dimensions = boundary.dimensions();
        double rz = a.coord.position().x(2);
        double dz1 = (0 + rz);
        double dz2 = (dimensions.x(2) - rz);
        return energy(dz1) + energy(dz2);
    }//end of energy

    private double energy(double r) {
        if (r < cutoff) {
            double rr = sigma / r;
            double r2 = rr * rr;
            double r6 = r2 * r2 * r2;
            return 4 * epsilon * r6 * (r6 - 1.0) + epsilon;
        }
        return 0;
    }

    private double gradient(double r) {
        if (r < cutoff) {
            double rr = sigma / r;
            double r2 = rr * rr;
            double r6 = r2 * r2 * r2;
            return -48 * epsilon * r6 * (r6 - 0.5);
        }
        return 0;
    }

    public Vector gradient(AtomSet atom) {
        Atom a = (Atom) atom;
        Vector dimensions = boundary.dimensions();
        double rz = a.coord.position().x(2);
        double dz1 = (dimensions.x(2) - rz);
        double gradz = gradient(rz) - gradient(dz1);
        gradient.setX(2, gradz);
        return gradient;
    }
    
    public double virial(AtomSet atoms) {
        return 0.0;
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
        cutoff = radius * Math.pow(2, 1. / 6.);
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
}