package etomica.virial;

import etomica.*;
import etomica.units.Dimension;

/**
 * @author kofke
 *
 * Meter to for perturbative calculation of a sum of cluster integrals.
 * 
 */
public class MeterVirial extends MeterGroup implements DatumSource {

	private final int N;
	private Cluster[] clusters;
	private Cluster refCluster;
	private double refIntegral;
	private double refTemperature, refBeta;
	private P2Cluster simulationPotential;
	private double temperature, beta;
	private PairSet pairSet;
	
	/**
	 * Constructor for MeterVirial.
	 * @param sim
	 */
	public MeterVirial(Simulation sim, PairSet pairSet, 
						double refTemperature, Cluster refCluster, double refIntegral, 
						Cluster[] clusters, P2Cluster simulationPotential) {
		super(sim, clusters.length+1);
		this.pairSet = pairSet;
		N = pairSet.nMolecules();
		this.clusters = clusters;
		this.refCluster = refCluster;
		this.refIntegral = refIntegral;
		this.simulationPotential = simulationPotential;
		setTemperature(Default.TEMPERATURE);
		setReferenceTemperature(refTemperature);
	}

	/**
	 * @see etomica.MeterGroup#updateValues()
	 */
	public void updateValues() {
		double pi = simulationPotential.pi();
		currentValues[0] = refCluster.value(pairSet, refBeta)/pi;
		for(int i=1; i<nMeters; i++) currentValues[i] = clusters[i-1].value(pairSet, beta)/pi;
	}
	
	public double value(DataSource.ValueType type) {
		MeterScalar[] m = allMeters();
		double sum = 0.0;
		for(int i=0; i<clusters.length; i++) {
			sum += clusters[i].weight()*m[i+1].average();
		}
		double denom = refCluster.weight()*m[0].average();
		return (denom != 0.0) ? sum*refIntegral/denom : Double.NaN;
	}

	/**
	 * @see etomica.MeterAbstract#getDimension()
	 */
	public Dimension getDimension() {
		return Dimension.NULL;
	}
	
	public void setPairs(PairSet pairs) {
		this.pairSet = pairs;
	}
	public PairSet getPairs() {
		return pairSet;
	}

	/**
	 * Sets the sigma.
	 * @param sigma The sigma to set
	 */
//	public void setSigma(double sigma) {
//		this.sigma = sigma;
//		sigma2 = sigma*sigma;
//		double b0 = 2.0*Math.PI/3.0*sigma2*sigma;
////		B4Ref = 0.3238*8*b0*b0*b0;
//	}
	
	/**
	 *  Overrides MCMoveAtom to prevent index-0 molecule from being displaced
	 */
	public static class MyMCMoveAtom extends MCMoveAtom {
		public MyMCMoveAtom(IntegratorMC integrator) {super(integrator);}
		
		public boolean doTrial() {
			if(phase.atomCount()==0) return false;
			atom = null;
//			double refPosition = phase.speciesMaster.atomList.getFirst().coord.position(0);
			while(atom == null || atom.node.index()==0) atom = phase.speciesMaster.atomList.getRandom();
			uOld = potential.calculate(phase, iteratorDirective.set(atom), energy.reset()).sum();
			atom.coord.displaceWithin(stepSize);
//			switch(atom.node.index()) {
//				case 2: atom.coord.position().setX(2,refPosition); //z = 0; fall through to case 1
//				case 1: atom.coord.position().setX(1,refPosition); //y = 0
//				default: 	
//			}
			uNew = Double.NaN;
			return true;
		}
	}

	/**
	 * Returns the temperature.
	 * @return double
	 */
	public double getTemperature() {
		return temperature;
	}

	/**
	 * Sets the temperature.
	 * @param temperature The temperature to set
	 */
	public void setTemperature(double temperature) {
		this.temperature = temperature;
		beta = 1.0/temperature;
	}

	/**
	 * Returns the refTemperature.
	 * @return double
	 */
	public double getReferenceTemperature() {
		return refTemperature;
	}

	/**
	 * Sets the refTemperature.
	 * @param refTemperature The refTemperature to set
	 */
	public void setReferenceTemperature(double refTemperature) {
		this.refTemperature = refTemperature;
		refBeta = 1.0/refTemperature;
	}

}
