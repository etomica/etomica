package etomica.virial;

import etomica.Phase;
import etomica.data.meter.MeterArray;
import etomica.units.Dimension;

/**
 * Measures value of clusters in a phase and returns the values
 * divided by the sampling bias of the integrator.
 */

public class MeterVirial extends MeterArray {

	protected final ClusterAbstract clusters[];
	double beta;
	protected final IntegratorClusterMC integrator;
    private double weightFactor;
	
	/**
	 * Constructor for MeterVirial.
	 */
	public MeterVirial(ClusterAbstract[] aClusters, IntegratorClusterMC aIntegrator, double temperature) {
		integrator = aIntegrator;
		clusters = aClusters;
		setTemperature(temperature);
		setNDataPerPhase(clusters.length);
        weightFactor = 1.0;
	}

    public void setPhase(Phase[] p) {
        if (p.length != 1) throw new IllegalArgumentException("MeterVirial can only handle 1 phase");
        super.setPhase(p);
        setSampleCluster(((PhaseCluster)p[0]).getSampleCluster());
    }

    /**
     * Sets the cluster used by the integrator to sample phase space.
     */
    public void setSampleCluster(ClusterWeight sampleCluster) {
        CoordinatePairSet cPairSet = ((PhaseCluster)phase[0]).getCPairSet();
        weightFactor = sampleCluster.value(cPairSet,beta) / integrator.getWeight();
    }
    
	public double[] getDataAsArray(Phase p) {
		CoordinatePairSet cPairSet = ((PhaseCluster)p).getCPairSet();
		double w = weightFactor * integrator.getWeight();
		for (int i=0; i<clusters.length; i++) {
			phaseData[i] = clusters[i].value(cPairSet,beta)/w;
		}
		return phaseData;
	}

	public Dimension getDimension() {
		return Dimension.FRACTION;
	}

	public double getTemperature() {
		return 1/beta;
	}

	public void setTemperature(double temperature) {
		beta = 1.0/temperature;
	}

}
