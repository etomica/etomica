/*
 * Created on Apr 9, 2004
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package etomica.modules.dcvgcmd;
import etomica.Atom;
import etomica.AtomSet;
import etomica.Default;
import etomica.EtomicaInfo;
import etomica.Space;
import etomica.potential.Potential1;
import etomica.potential.PotentialSoft;
import etomica.space.Vector;
/**
 * @author Owner
 *
 * To change the template for this generated type comment go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */

/*!!!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!
 * !!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!
 * Radius is actually sigma.  Be forewarned!!!!
 * @author ecc4
 *
 * To change the template for this generated type comment go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */

public class P1LJWCAWall extends Potential1 implements PotentialSoft{
	
	private final int D;
		private final Vector gradient;
		private double radius, radius2;
		private double cutoff;
		private double cutoff2;
	
		public P1LJWCAWall(Space parent) {
			super(parent);
			D = space.D();
			gradient = space.makeVector();
			setRadius(0.5*Default.ATOM_SIZE);
		}
    
		public static EtomicaInfo getEtomicaInfo() {
			EtomicaInfo info = new EtomicaInfo("WCA LJ Potential in the Z-Coordinate");
			return info;
		}
		
		public double energy(AtomSet atom) {
			Atom a = (Atom)atom;
			Vector dimensions = a.node.parentPhase().boundary().dimensions();
			double rz = a.coord.position().x(2);
			double dz1 = (0 + rz);
			double dz2 = (dimensions.x(2) - rz);
			return energy(dz1) + energy(dz2);		
		}//end of energy
	
		private double energy(double r) {
			double rr = radius/r;
			double r2 = rr*rr;
			double r6 = r2*r2*r2;
			if(r*r < cutoff2) {
                return 4*Default.POTENTIAL_WELL*r6*(r6 - 1.0)+Default.POTENTIAL_WELL;
			}
            return 0;
		}
	
		private double gradient(double r) {
			double rr = radius/r;
			double r2 = rr*rr;
			double r6 = r2*r2*r2;
			if(r*r < cutoff2) {
			    return -48*Default.POTENTIAL_WELL*r6*(r6 - 0.5);
			}
			return 0;
		}
	
		public Vector gradient(AtomSet atom) {
			Atom a = (Atom)atom;
			Vector dimensions = a.node.parentPhase().boundary().dimensions();
			double rz = a.coord.position().x(2);
			double dz1 = (dimensions.x(2) - rz);
			double gradz = gradient(rz) - gradient(dz1);
			gradient.setX(2,gradz);
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
            cutoff = radius*Math.pow(2,1./6.);
            cutoff2 = cutoff*cutoff;
		}

}	


