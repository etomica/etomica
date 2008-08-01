package etomica.virial;

import etomica.data.Data;
import etomica.data.DataSource;
import etomica.data.DataSourceAtomDistance;
import etomica.data.DataTag;
import etomica.data.IDataInfo;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.units.Null;

/**
 * Measures value of clusters in a box and returns the values
 * divided by the sampling bias from the sampling cluster.
 */
public class MeterSamplingWeight implements DataSource, java.io.Serializable {

    /**
	 * Constructor for MeterVirial.
	 */
	public MeterSamplingWeight(DataSourceAtomDistance dataDistance) {
		this.dataDistance = dataDistance; 
        data = new DataDoubleArray(2);
        dataInfo = new DataInfoDoubleArray("Cluster Value",Null.DIMENSION, new int[]{2});
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
        double pi = box.getSampleCluster().value(box);
        data.getData()[0] = dataDistance.getDataAsScalar();
        data.getData()[1] = pi;
        return data;
    }
    
    public BoxCluster getBox() {
        return box;
    }
    
    public void setBox(BoxCluster newBox) {
        box = newBox;
    }

    protected final DataSourceAtomDistance dataDistance;
	private final DataDoubleArray data;
	private final IDataInfo dataInfo;
    private final DataTag tag;
    private BoxCluster box;
    private static final long serialVersionUID = 1L;
}
