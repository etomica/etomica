/*
 * Created on Apr 9, 2004
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package etomica.modules.dcvgcmd;

import etomica.EtomicaInfo;
import etomica.atom.Atom;
import etomica.atom.AtomSet;
import etomica.potential.Potential1;
import etomica.potential.PotentialSoft;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * @author Owner
 * 
 * To change the template for this generated type comment go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */

public class P1WCAPorousWall extends Potential1 implements PotentialSoft {

    private final Vector gradient;
    private double sigma, sigma2;
    private double epsilon;
    private double cutoff, cutoff2;
    private double poreRadius, poreRadius2;
    private Vector[] poreCenters;
    private double z;

    public P1WCAPorousWall(Simulation sim) {
        this(sim.space, sim.getDefaults().atomSize, sim.getDefaults().potentialWell);
    }

    public P1WCAPorousWall(Space space, double sigma, double epsilon) {
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
        Vector r = a.coord.position();
        double rz = r.x(2);
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
    
    public double virial(AtomSet atoms) {
        return 0.0;
    }
    
    private boolean inPore(Vector r) {
        for(int i=0; i<poreCenters.length; i++) {
            double dx = r.x(0) - poreCenters[i].x(0);
            double dy = r.x(1) - poreCenters[i].x(1);
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

    public Vector gradient(AtomSet atom) {
        Atom a = (Atom) atom;
        Vector r = a.coord.position();
        double rz = r.x(2);
        double dz2 = (z - rz);
        dz2 *= dz2;
        double gradz = 0.0;
        if(dz2 < cutoff2 && !inPore(r)) {
            gradz = gradient(dz2);
            if(z > rz) gradz = -gradz;
//            System.out.println("rz, grad: "+rz + "  "+gradz);
        }
        gradient.setX(2, gradz);
        return gradient;
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