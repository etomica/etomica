package etomica.research.nonequilwork;

import etomica.*;

/**
 * @author kofke
 *
 * Overlap-sampling potential in which the two systems differ in the presence of
 * a single molecule.  Applicable only in simple case of a single species (one
 * pair potential).  Cannot compute quantities (e.g., virial) other than the
 * energy.
 */
public class PotentialOSInsert extends PotentialN {

	public PotentialOSInsert(SimulationElement parent) {
		super(parent);
		atomPair = new AtomPair(simulation().space);
	}

	/**
	 * @see etomica.Potential#calculate(etomica.AtomSet, etomica.IteratorDirective, etomica.PotentialCalculation)
	 */
	//basis should be the SpeciesAgent corresponding to the set Species, and it is just passed on through to the pair potential
	public void calculate(AtomSet basis, IteratorDirective id, PotentialCalculation pc) {
		if(!enabled) return;
		if( !(pc instanceof PotentialCalculationEnergySum) ) throw new IllegalArgumentException("PotentialOSInsert suitable only for energy calculation");

		Atom atom1 = id.atom1();
		
		potential.calculate(basis, testDirective, testEnergySum);
		uTest = testEnergySum.sum();
		
		if(id.atomCount() == 0) {
			potential.calculate(basis, id, energySum);//total (N+1)-interaction energy
			u0 = energySum.sum() - uTest;//N-interaction energy, without test molecule
		} else if(atom1 == testMolecule) {
			u0 = 0.0;//only test molecule energy involved
		} else {
			potential.calculate(basis, id, energySum);//atom1 energy sum
			double u0T = potential.energy(atomPair.reset(testMolecule, atom1));//atom1-testMolecule pair energy
			u0 = energySum.sum() - u0T;//atom1 energy without interation with testMolecule
		}
		
 		updateW();
		pc.set(this).actionPerformed((AtomSet)null);
	}
	
	/**
	 * Argument is ignored.  Energy is given by terms computed in last call to
	 * calculate method.
	 * @see etomica.PotentialN#energy(AtomSet)
	 */
	public double energy(AtomSet atomN) {
		return u0 + W;
	}

	/**
	 * @see etomica.Potential#setIterator(etomica.AtomSetIterator)
	 */
	public void setIterator(AtomSetIterator iterator) {
		throw new etomica.exception.MethodNotImplementedException();
	}

	/**
	 * @see etomica.Potential#getIterator()
	 */
	public AtomSetIterator getIterator() {
		throw new etomica.exception.MethodNotImplementedException();
	}
	
	public void setSpecies(Species species) {
		super.setSpecies(new Species[] {species});
		potential.setIterator(new Api1A(new ApiIntragroup1A(simulation()),
								new ApiIntragroupAA(simulation())));
	}

	/**
	 * Returns the potential.
	 * @return Potential2
	 */
	public Potential2 getPotential() {
		return potential;
	}

	/**
	 * Sets the potential.
	 * @param potential The potential to set
	 */
	public void setPotential(Potential2 potential) {
		this.potential = potential;
	}
	
	public void setGamma(double gamma){ 
		this.gamma = gamma; 
		if(gamma == 0.0) Li = Double.NEGATIVE_INFINITY;
		else Li = Math.log(gamma);
	}
	public void setLi(double Li){ 
		this.Li = Li; 
		if(Double.isInfinite(Li)) gamma = 0.0;
		else gamma = Math.exp(Li);
	}
	public double getGamma(){ return gamma;}
	public double getW() {
		updateW();
//		System.out.println("Utest "+uTest+" gamma "+gamma+" Li "+Li);
		return W;
	}        
	private void updateW() {
		if(Double.isInfinite(Li)) W = 0;
		else if(Double.isInfinite(uTest)) W = Double.POSITIVE_INFINITY;
		else W = Math.log((1-gamma)+Math.exp(beta*(uTest-C)+Li))/beta;
	}
	public double getUTest() {return uTest;}
	public double getCValue() {return C;}
	public void setCValue(double C) {this.C = C;}


	private final IteratorDirective testDirective = new IteratorDirective(IteratorDirective.BOTH);
	private final PotentialCalculationEnergySum energySum = new PotentialCalculationEnergySum();
	private final PotentialCalculationEnergySum testEnergySum = new PotentialCalculationEnergySum();
	private Atom testMolecule;
	private double gamma;
	private double Li = Double.NEGATIVE_INFINITY;//ln(gamma) = Li.
	private double beta;
	private double W, uTest, u0;
	private double C;
	private final AtomPair atomPair;
	private Potential2 potential;


	/**
	 * Returns the testMolecule.
	 * @return Atom
	 */
	public Atom getTestMolecule() {
		return testMolecule;
	}

	/**
	 * Sets the testMolecule.
	 * @param testMolecule The testMolecule to set
	 */
	public void setTestMolecule(Atom testMolecule) {
		this.testMolecule = testMolecule;
		testDirective.set(testMolecule);
	}

}
