/*
 * Created on Jan 30, 2004

 */
package etomica.chem.models.water;

import etomica.chem.Model;
import etomica.chem.electrostatics.Charge;
import etomica.chem.elements.Hydrogen;
import etomica.chem.elements.Oxygen;
import etomica.chem.models.*;
import etomica.*;

/**
 * Generic model for a 3-site representation of water, with one Lennard-
 * Jones oxygen site and two point-charge hydrogen sites.  Subclasses differ in
 * the LJ parameters, point charges, and geometry of the model.
 */
public class Abstract3Site extends ModelMolecular {

	/**
	 * @param sigma Lennard-Jones sigma parameter for oxygen-oxygen interaction
	 * @param epsilon Lennard-Jones epsilon parameter for oxygen-oxygen
	 * interaction
	 * @param qH partial charge on the hydrogen; oxygen charge will be -2qH
	 * @param rOH oxygen-hydrogen bond length
	 * @param theta H-O-H angle
	 */
	public Abstract3Site(double sigma, double epsilon, double qH, double rOH, double theta) {
		super(new ConfigurationWater(rOH, theta), makeModels(sigma, epsilon, qH), new int[] {1, 2});
	}

	private static Model[] makeModels(double sigma, double epsilon, double qH) {
		ModelAtomic oxygen = new LennardJones(Oxygen.INSTANCE, sigma, epsilon, new Charge(-2*qH));
		ModelAtomic hydrogen = new Point(Hydrogen.INSTANCE, new Charge(qH));
		return new Model[] {oxygen, hydrogen};
	}
	
	public Potential makePotential(Space space) {
		PotentialWW potential = new PotentialWW(space, truncation);
		potential.setParameters((LennardJones)models[0], (Charge)((ModelAtomic)models[1]).getElectrostatic());
		return potential;
	}
	
	public static void main(String[] args) {
		etomica.graphics.SimulationGraphic sim = new etomica.graphics.SimulationGraphic(Space3D.INSTANCE);
		Model waterModel = new Abstract3Site(3.0, 100.0, 0.5, 1.0, Math.PI*109./180.);
		Species species = new Species(Space3D.INSTANCE, waterModel);
		Phase phase = new Phase(Space3D.INSTANCE);
		etomica.graphics.DisplayPhase display = new etomica.graphics.DisplayPhase(sim);
		sim.elementCoordinator.go();
		sim.makeAndDisplayFrame();
	}
	
	public static class PotentialWW extends Potential2 implements Potential2.Soft {

		public PotentialWW(Space space, PotentialTruncation potentialTruncation) {
			super(space,potentialTruncation);
			dr = (Space3D.Vector)space.makeVector();
			shift = (Space3D.Vector)space.makeVector();
		}   
		public double energy(Atom[] pair){
			double sum = 0.0;
			double r2 = 0.0;
			
			Atom[] atomArray1 = ((AtomTreeNodeGroupArray)pair[0].node).childAtomArray();
			Atom[] atomArray2 = ((AtomTreeNodeGroupArray)pair[1].node).childAtomArray();
		
			//compute O-O distance to consider truncation	
			Space3D.Vector O1r = (Space3D.Vector)atomArray1[0].coord.position();
			Space3D.Vector O2r = (Space3D.Vector)atomArray2[0].coord.position();

			dr.Ev1Mv2(O1r, O2r);
			boundary.nearestImage(dr, shift);
			r2 = dr.squared();

			if(potentialTruncation.isZero(r2)) return 0.0;

			if(r2<1.6) return Double.POSITIVE_INFINITY;
		
			
			sum += chargeOO/Math.sqrt(r2);
			double s2 = sigma2/r2;
			double s6 = s2*s2*s2;
			sum += epsilon4*s6*(s6 - 1.0);
		
			Space3D.Vector H11r = (Space3D.Vector)atomArray1[1].coord.position();
			Space3D.Vector H12r = (Space3D.Vector)atomArray1[2].coord.position();
			Space3D.Vector H21r = (Space3D.Vector)atomArray2[1].coord.position();
			Space3D.Vector H22r = (Space3D.Vector)atomArray2[2].coord.position();
        		
			final boolean zeroShift = shift.isZero();
					
			r2 = (zeroShift) ? O1r.Mv1Squared(H21r) : O1r.Mv1Pv2Squared(H21r,shift);
			if(r2<1.6) return Double.POSITIVE_INFINITY;
			sum += chargeOH/Math.sqrt(r2);
		
			r2 = (zeroShift) ? O1r.Mv1Squared(H22r) : O1r.Mv1Pv2Squared(H22r,shift);
			if(r2<1.6) return Double.POSITIVE_INFINITY;
			sum += chargeOH/Math.sqrt(r2);

			r2 = (zeroShift) ? H11r.Mv1Squared(O2r) : H11r.Mv1Pv2Squared(O2r,shift);
			if(r2<1.6) return Double.POSITIVE_INFINITY;
			sum += chargeOH/Math.sqrt(r2);

			r2 = (zeroShift) ? H11r.Mv1Squared(H21r) : H11r.Mv1Pv2Squared(H21r,shift);
			if(r2<1.6) return Double.POSITIVE_INFINITY;
			sum += chargeHH/Math.sqrt(r2);

			r2 = (zeroShift) ? H11r.Mv1Squared(H22r) : H11r.Mv1Pv2Squared(H22r,shift);
			if(r2<1.6) return Double.POSITIVE_INFINITY;
			sum += chargeHH/Math.sqrt(r2);

			r2 = (zeroShift) ? H12r.Mv1Squared(O2r) : H12r.Mv1Pv2Squared(O2r,shift);
			if(r2<1.6) return Double.POSITIVE_INFINITY;
			sum += chargeOH/Math.sqrt(r2);

			r2 = (zeroShift) ? H12r.Mv1Squared(H21r) : H12r.Mv1Pv2Squared(H21r,shift);
			if(r2<1.6) return Double.POSITIVE_INFINITY;
			sum += chargeHH/Math.sqrt(r2);

			r2 = (zeroShift) ? H12r.Mv1Squared(H22r) : H12r.Mv1Pv2Squared(H22r,shift);
			if(r2<1.6) return Double.POSITIVE_INFINITY;
			sum += chargeHH/Math.sqrt(r2);

			return sum;																					        
		}//end of energy
    
		public Space.Vector gradient(Atom[] pair){
			throw new etomica.exception.MethodNotImplementedException();
		}
		public double hyperVirial(Atom[] pair){
			throw new etomica.exception.MethodNotImplementedException();
		}
		public double integral(double rC){
			throw new etomica.exception.MethodNotImplementedException();
		}
		public double virial(Atom[] pair){
			throw new etomica.exception.MethodNotImplementedException();
		}
    
		public double getSigma() {return sigma;}    
		public double getEpsilon() {return epsilon;}
    
		private final void setParameters(LennardJones lj, Charge qH) {
			sigma = lj.getSigma();
			sigma2 = sigma*sigma;
			epsilon = lj.getEpsilon();
			epsilon4 = 4*epsilon;
			chargeH = qH.getQ();
			chargeO = -2.0*chargeH;
			chargeOO = chargeO * chargeO;
			chargeOH = chargeO * chargeH;
			chargeHH = chargeH * chargeH;
		}
    
		public double sigma , sigma2;
		public double epsilon, epsilon4;
		private Space3D.Boundary boundary;
		private double chargeH;
		private double chargeO;
		private double chargeOO, chargeOH, chargeHH;
		private Space3D.Vector dr, shift;
		/**
		 * Returns the boundary.
		 * @return Space3D.Boundary
		 */
		public Space3D.Boundary getBoundary() {
			return boundary;
		}

		/**
		 * Sets the boundary.
		 * @param boundary The boundary to set
		 */
		public void setBoundary(Space3D.Boundary boundary) {
			this.boundary = boundary;
		}
	}//end of PotentialWW
	
	public static class ConfigurationWater extends etomica.Configuration {

		private double bondLengthOH = 1.0;
		private double angleHOH = 109.5*Math.PI/180.;

		public ConfigurationWater(Space space) {
			super(space);
		}
    
		public ConfigurationWater(double rOH, double theta) {
			super(new Space3D());
			bondLengthOH = rOH;
			angleHOH = theta;
		}
    
		public void initializePositions(AtomIterator[] iterators){
			if(iterators == null || iterators.length == 0) return;
        
			AtomIterator iterator = iterators[0];
			double x = 0.0;
			double y = 0.0;
        
			iterator.reset();
        
			Atom o = iterator.nextAtom();
			o.coord.position().E(new double[] {x, y, 0.0});
               
			Atom h1 = iterator.nextAtom();
			h1.coord.position().E(new double[] {x+bondLengthOH, y, 0.0});
                
			Atom h2 = iterator.nextAtom();
			h2.coord.position().E(new double[] {x+bondLengthOH*Math.cos(angleHOH), y+bondLengthOH*Math.sin(angleHOH), 0.0});

		}//end of initializePositions
	}//end of Configuration
}
