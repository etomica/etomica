package etomica.virial;

import etomica.DataSource;
import etomica.Simulation;
import etomica.data.meter.MeterGroup;
import etomica.data.meter.MeterScalar;
import etomica.units.Dimension;

/**
 * @author kofke
 *
 * Computes cluster integral using sampling cluster-potential that is the same
 * as the target cluster.  
 * Gamma = Gamma0 <abs(gamma)>_gamma / <gamma0/gamma>_gamma
 */

/* History
 * 12/16/03 (DAK) new, from MeterVirial
 */
public class MeterVirialAbs extends MeterGroup implements DatumSource {

	private final int N;
	private ClusterAbstract refCluster;
	private double refIntegral;
	private double refTemperature, refBeta;
	private P0Cluster simulationPotential;
	
	/**
	 * Constructor for MeterVirial.
	 * @param sim
	 */
	public MeterVirialAbs(Simulation sim, 
						double refTemperature, ClusterAbstract refCluster, double refIntegral, 
						P0Cluster simulationPotential) {
		super(sim, 2);
		N = refCluster.pointCount();
		this.refCluster = refCluster;
		this.refIntegral = refIntegral;
		this.simulationPotential = simulationPotential;
		setReferenceTemperature(refTemperature);
		allMeters()[0].setLabel("Ref perturb");
		for(int i=1; i<nMeters; i++) {
			allMeters()[i].setLabel("Cluster "+i);
		}
	}

	/**
	 * @see etomica.data.meter.MeterGroup#updateValues()
	 */
	public void updateValues() {
		double pi = simulationPotential.pi((PhaseCluster)phase);
		PairSet pairSet = ((PhaseCluster)phase).getPairSet().resetPairs();//resetPairs not needed if can be sure it is done in potential.pi method call
		currentValues[0] = refCluster.value(pairSet, refBeta)/pi;
		currentValues[1] = simulationPotential.mostRecentIsPositive() ? +1 : -1;
	}
	
	public double value(DataSource.ValueType type) {
		MeterScalar[] m = allMeters();
		double sum = simulationPotential.getCluster().weight()*m[1].average();
		double denom = refCluster.weight()*m[0].average();
		return (denom != 0.0) ? sum*refIntegral/denom : Double.NaN;
	}

	/**
	 * @see etomica.MeterAbstract#getDimension()
	 */
	public Dimension getDimension() {
		return Dimension.NULL;
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
