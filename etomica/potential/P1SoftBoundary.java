package etomica.potential;

import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.atom.Atom;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomSet;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.Dimension;
import etomica.units.Length;

/**
 * @author David Kofke
 *
 * Inverse-power potential between an atom and all four boundaries of the phase.  Potential
 * is of the form <tt>u(r) = (R/r)^12</tt>, where <tt>R</tt> is the repulsion radius.  This
 * term is summed over all four boundaries.
 */
public class P1SoftBoundary extends Potential1 implements PotentialSoft, EtomicaElement {

    private static final long serialVersionUID = 1L;
	private final Vector[] gradient;
	private double radius;
	
    public P1SoftBoundary(Simulation sim) {
        this(sim.getSpace(), 0.5*sim.getDefaults().atomSize);
    }
	public P1SoftBoundary(Space space, double radius) {
		super(space);
        gradient = new Vector[1];
		gradient[0] = space.makeVector();
		setRadius(radius);
	}
    
	public static EtomicaInfo getEtomicaInfo() {
		EtomicaInfo info = new EtomicaInfo("PotentialSoft repulsive potential at the phase boundaries");
		return info;
	}
    
	public double energy(AtomSet a) {
		Vector dimensions = ((Atom)a).getNode().parentPhase().getBoundary().getDimensions();
		double rx = ((AtomLeaf)a).getCoord().getPosition().x(0);
		double ry = ((AtomLeaf)a).getCoord().getPosition().x(1);
		double dx1 = (dimensions.x(0) - rx);
		double dy1 = (dimensions.x(1) - ry);
		return energy(rx) + energy(ry) + energy(dx1) + energy(dy1);		
	}//end of energy
	
	private double energy(double r) {
		r /= radius;
		double r2 = 1./(r*r);
		double r6 = r2*r2*r2;
		return r6*r6;
	}
	
	private double gradient(double r) {
		double rr = radius/r;
		double r2 = rr*rr;
		double r6 = r2*r2*r2;
		return -12*r6*r6/r;
	}
	
	public Vector[] gradient(AtomSet a) {
		Vector dimensions = boundary.getDimensions();
		double rx = ((AtomLeaf)a).getCoord().getPosition().x(0);
		double ry = ((AtomLeaf)a).getCoord().getPosition().x(1);
		double dx1 = (dimensions.x(0) - rx);
		double dy1 = (dimensions.x(1) - ry);
		double gradx = gradient(rx) - gradient(dx1);
		double grady = gradient(ry) - gradient(dy1);
		gradient[0].setX(0,gradx);
		gradient[0].setX(1,grady);
		return gradient;
	}
	
	public double virial(AtomSet atoms) {
	    return 0.0;
    }
    
	/**
	 * Returns the radius.
	 * @return double
	 */
	public double getRadius() {
		return radius;
	}

	/**
	 * Sets the radius.
	 * @param radius The radius to set
	 */
	public void setRadius(double radius) {
		this.radius = radius;
	}
    
    public Dimension getRadiusDimension() {
        return Length.DIMENSION;
    }
    
}
