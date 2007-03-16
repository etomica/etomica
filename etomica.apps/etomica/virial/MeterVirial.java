package etomica.virial;

import etomica.data.Data;
import etomica.data.DataSource;
import etomica.data.DataTag;
import etomica.data.IDataInfo;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.units.Null;

/**
 * Measures value of clusters in a phase and returns the values
 * divided by the sampling bias from the sampling cluster.
 */
public class MeterVirial implements DataSource, java.io.Serializable {

    /**
	 * Constructor for MeterVirial.
	 */
	public MeterVirial(ClusterAbstract[] aClusters) {
		clusters = aClusters;
        data = new DataDoubleArray(clusters.length);
        dataInfo = new DataInfoDoubleArray("Cluster Value",Null.DIMENSION, new int[]{clusters.length});
        tag = new DataTag();
        dataInfo.addTag(tag);
	}

	public IDataInfo getDataInfo() {
        return dataInfo;
    }
    
    public DataTag getTag() {
        return tag;
    }
    
    public Data getData() {
        CoordinatePairSet cPairSet = phase.getCPairSet();
        AtomPairSet aPairSet = phase.getAPairSet();
        double pi = phase.getSampleCluster().value(cPairSet,aPairSet);
        double x[] = data.getData();
        for (int i=0; i<clusters.length; i++) {
            x[i] = clusters[i].value(cPairSet,aPairSet)/pi;
        }
        return data;
    }
    
    public PhaseCluster getPhase() {
        return phase;
    }
    
    public void setPhase(PhaseCluster newPhase) {
        phase = newPhase;
    }

    protected final ClusterAbstract clusters[];
	private final DataDoubleArray data;
	private final IDataInfo dataInfo;
    private final DataTag tag;
    private PhaseCluster phase;
    private static final long serialVersionUID = 1L;
}
