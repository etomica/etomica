package etomica.virial;

import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.units.Null;

/**
 * Measures value of clusters in a box and returns the values
 * divided by the sampling bias from the sampling cluster.
 */
public class MeterVirial implements IEtomicaDataSource, java.io.Serializable {

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

	public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }
    
    public DataTag getTag() {
        return tag;
    }
    
    public IData getData() {
        double pi = box.getSampleCluster().value(box);
        double x[] = data.getData();
        for (int i=0; i<clusters.length; i++) {
            x[i] = clusters[i].value(box)/pi;
            if (Double.isNaN(x[i])) throw new RuntimeException("oops "+clusters[i].value(box)+" "+x[i]+" "+pi);
        }
        return data;
    }
    
    public ClusterAbstract[] getClusters() {
        return clusters;
    }
    
    public BoxCluster getBox() {
        return box;
    }
    
    public void setBox(BoxCluster newBox) {
        box = newBox;
    }

    protected final ClusterAbstract clusters[];
	private final DataDoubleArray data;
	private final IEtomicaDataInfo dataInfo;
    private final DataTag tag;
    private BoxCluster box;
    private static final long serialVersionUID = 1L;
}
