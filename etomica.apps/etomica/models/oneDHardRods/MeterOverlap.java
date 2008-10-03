package etomica.models.oneDHardRods;

import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.DataSource;
import etomica.data.DataSourceScalar;
import etomica.data.DataTag;
import etomica.data.IDataInfo;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.units.Dimension;


public class MeterOverlap implements DataSource {
    DataTag tag;
    DataInfo dataInfo;
    DataSourceScalar dataSourceA, dataSourceB;
    double temperature;
    DataDoubleArray dda;
    
    MeterOverlap(String label, Dimension dimension, DataSourceScalar dataSourceA,
            DataSourceScalar dataSourceB, double temperature){
        dataInfo = new DataInfoDouble(label, dimension);
        tag = new DataTag();
        dataInfo.addTag(tag);
        
        this.dataSourceA = dataSourceA;
        this.dataSourceB = dataSourceB;
        this.temperature = temperature;
    }
    
    public DataDoubleArray getData(){
        double[] eAeB = dda.getData();
        
        eAeB[1] = Math.exp(-dataSourceB.getDataAsScalar()/temperature);
        eAeB[0] = Math.exp(-dataSourceA.getDataAsScalar()/temperature);
        
        return dda;
    }
    
    public IDataInfo getDataInfo(){
        return dataInfo;
    }
    public DataTag getTag() {
        return tag;
    }
}
