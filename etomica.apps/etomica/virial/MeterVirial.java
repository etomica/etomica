package etomica.virial;

import etomica.Phase;
import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.DataSource;
import etomica.data.Meter;
import etomica.data.types.DataDoubleArray;
import etomica.units.Dimension;
import etomica.utility.NameMaker;

/**
 * Measures value of clusters in a phase and returns the values
 * divided by the sampling bias of the integrator.
 */

public class MeterVirial implements DataSource, Meter, java.io.Serializable {

	protected final ClusterAbstract clusters[];
	protected final IntegratorClusterMC integrator;
    private double weightFactor;
	
	/**
	 * Constructor for MeterVirial.
	 */
	public MeterVirial(ClusterAbstract[] aClusters, IntegratorClusterMC aIntegrator) {
        setName(NameMaker.makeName(this.getClass()));
		integrator = aIntegrator;
		clusters = aClusters;
        data = new DataDoubleArray("Cluster Value",Dimension.NULL,clusters.length);
        weightFactor = 1.0;
	}

	public DataInfo getDataInfo() {
        return data.getDataInfo();
    }
    
    /**
     * Sets the cluster used by the integrator to sample phase space.
     */
    //XXX this is really fragile.  If the integrator weight is != 1, resetting the 
    //integrator (which sets the weight to 1) after calling this method
    //breaks the class
    public void setSampleCluster(ClusterWeight sampleCluster) {
        CoordinatePairSet cPairSet = phase.getCPairSet();
        AtomPairSet aPairSet = phase.getAPairSet();
        weightFactor = sampleCluster.value(cPairSet,aPairSet) / integrator.getWeight();
    }
    
	public Data getData() {
		CoordinatePairSet cPairSet = phase.getCPairSet();
        AtomPairSet aPairSet = phase.getAPairSet();
		double w = weightFactor * integrator.getWeight();
        double[] x = data.getData();
		for (int i=0; i<clusters.length; i++) {
			x[i] = clusters[i].value(cPairSet,aPairSet)/w;
//            System.out.println("in meter, v="+x[i]*w+" w="+w+" "+integrator);
		}
		return data;
	}

	public Dimension getDimension() {
		return Dimension.FRACTION;
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
        this.phase = (PhaseCluster)phase;
        setSampleCluster(((PhaseCluster)phase).getSampleCluster());
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
