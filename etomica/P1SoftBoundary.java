package etomica;

/**
 * @author David Kofke
 *
 * Inverse-power potential between an atom and all four boundaries of the phase.
 */
public class P1SoftBoundary extends Potential1 implements PotentialSoft, EtomicaElement {

	private final int D;
	private final Space.Vector gradient;
	private double radius, radius2;
	private double cutoff = Double.MAX_VALUE;
	private Atom atom;
	
	public P1SoftBoundary() {
		this(Simulation.instance.hamiltonian.potential);
	}
    
	public P1SoftBoundary(SimulationElement parent) {
		super(parent);
		D = space.D();
		gradient = space.makeVector();
		setRadius(0.5*Default.ATOM_SIZE);
	}
    
	public static EtomicaInfo getEtomicaInfo() {
		EtomicaInfo info = new EtomicaInfo("PotentialSoft repulsive potential at the phase boundaries");
		return info;
	}
    
	public double energy(Atom[] a) {
		atom = a[0];
		Space.Vector dimensions = atom.node.parentPhase().dimensions();
		double rx = atom.coord.position(0);
		double ry = atom.coord.position(1);
		double dx1 = (dimensions.x(0) - rx);
		double dy1 = (dimensions.x(1) - ry);
		return energy(rx) + energy(ry) + energy(dx1) + energy(dy1);		
	}//end of energy
	
	private double energy(double r) {
		r /= radius;
		double r2 = r*r;
		double r6 = r2*r2*r2;
		return r6*r6;
	}
	
	private double gradient(double r) {
		double rr = radius/r;
		double r2 = rr*rr;
		double r6 = r2*r2*r2;
		return -12*r6*r6/r;
	}
	
	public Space.Vector gradient(Atom[] a) {
		atom = a[0];
		Space.Vector dimensions = boundary.dimensions();
		double rx = atom.coord.position(0);
		double ry = atom.coord.position(1);
		double dx1 = (dimensions.x(0) - rx);
		double dy1 = (dimensions.x(1) - ry);
		double gradx = gradient(rx) - gradient(dx1);
		double grady = gradient(ry) - gradient(dy1);
		gradient.setX(0,gradx);
		gradient.setX(1,grady);
		return gradient;
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
		radius2 = radius*radius;
	}

}
