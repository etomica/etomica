package etomica.virial;

import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.DataSource;
import etomica.data.meter.Meter;
import etomica.data.types.DataDoubleArray;
import etomica.integrator.IntegratorMC;
import etomica.phase.Phase;
import etomica.units.Dimension;
import etomica.units.Fraction;
import etomica.units.Null;
import etomica.util.NameMaker;

/**
 * Measures value of clusters in a phase and returns the values
 * divided by the sampling bias of the integrator.
 */

public class MeterVirial implements DataSource, Meter, java.io.Serializable {

	protected final ClusterAbstract clusters[];
	protected final IntegratorMC integrator;
	
	/**
	 * Constructor for MeterVirial.
	 */
	public MeterVirial(ClusterAbstract[] aClusters, IntegratorMC aIntegrator) {
        setName(NameMaker.makeName(this.getClass()));
		integrator = aIntegrator;
		clusters = aClusters;
        data = new DataDoubleArray("Cluster Value",Null.DIMENSION,clusters.length);
	}

	public DataInfo getDataInfo() {
        return data.getDataInfo();
    }
    
	public Data getData() {
		CoordinatePairSet cPairSet = phase.getCPairSet();
        AtomPairSet aPairSet = phase.getAPairSet();
        double w = sampleCluster.value(cPairSet,aPairSet);
        double[] x = data.getData();
		for (int i=0; i<clusters.length; i++) {
			x[i] = clusters[i].value(cPairSet,aPairSet)/w;
//            System.out.println("in meter "+getName()+" "+i+" "+clusters[i]+", v="+x[i]*w+" w="+w+" "+integrator);
		}
		return data;
	}

	public Dimension getDimension() {
		return Fraction.DIMENSION;
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
        sampleCluster = ((PhaseCluster)phase).getSampleCluster();
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
    private ClusterWeight sampleCluster;
}
