package etomica.virial;

import etomica.Data;
import etomica.DataInfo;
import etomica.DataSource;
import etomica.Meter;
import etomica.Phase;
import etomica.data.types.DataDoubleArray;
import etomica.units.Dimension;
import etomica.utility.NameMaker;

/**
 * Measures value of clusters in a phase and returns the values
 * divided by the sampling bias of the integrator.
 */

public class MeterVirial implements DataSource, Meter, java.io.Serializable {

	protected final ClusterAbstract clusters[];
	double beta;
	protected final IntegratorClusterMC integrator;
    private double weightFactor;
	
	/**
	 * Constructor for MeterVirial.
	 */
	public MeterVirial(ClusterAbstract[] aClusters, IntegratorClusterMC aIntegrator, double temperature) {
        setName(NameMaker.makeName(this.getClass()));
		integrator = aIntegrator;
		clusters = aClusters;
		setTemperature(temperature);
        data = new DataDoubleArray(new DataInfo("Cluster Value",Dimension.NULL));
		data.setLength(clusters.length);
        weightFactor = 1.0;
	}

	public DataInfo getDataInfo() {
        return data.getDataInfo();
    }
    
    /**
     * Sets the cluster used by the integrator to sample phase space.
     */
    public void setSampleCluster(ClusterWeight sampleCluster) {
        CoordinatePairSet cPairSet = phase.getCPairSet();
        AtomPairSet aPairSet = phase.getAPairSet();
        weightFactor = sampleCluster.value(cPairSet,aPairSet,beta) / integrator.getWeight();
    }
    
	public Data getData() {
		CoordinatePairSet cPairSet = phase.getCPairSet();
        AtomPairSet aPairSet = phase.getAPairSet();
		double w = weightFactor * integrator.getWeight();
        double[] x = data.getData();
		for (int i=0; i<clusters.length; i++) {
			x[i] = clusters[i].value(cPairSet,aPairSet,beta)/w;
		}
		return data;
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

    /**
     * @return Returns the phase.
     */
    public Phase getPhase() {
        return phase;
    }
    /**
     * @param phase The phase to set.
     */
    public void setPhase(Phase phase) {
        setSampleCluster(((PhaseCluster)phase).getSampleCluster());
        this.phase = (PhaseCluster)phase;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }
    private String name;
    private PhaseCluster phase;
	private final DataDoubleArray data;
}
