package etomica.models.oneDHardRods;

import etomica.data.DataInfo;
import etomica.data.DataSource;
import etomica.data.DataSourceScalar;
import etomica.data.DataTag;
import etomica.data.IDataInfo;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.units.Dimension;


public class MeterOverlap implements DataSource {
    DataTag tag;
    DataInfo dataInfo;
    DataSourceScalar dataSourceA, dataSourceB;
    double temperature;
    DataDoubleArray dda;
    
    /**
     * Put the system you are measuring in as the first DataSourceScalar
     * @param label
     * @param dimension
     * @param dataSourceA - denominator
     * @param dataSourceB - numerator
     * @param temperature
     */
    MeterOverlap(String label, Dimension dimension, DataSourceScalar dataSourceA,
            DataSourceScalar dataSourceB, double temperature){
        dataInfo = new DataInfoDoubleArray(label, dimension, new int[]{2});
        tag = new DataTag();
        dataInfo.addTag(tag);
        
        this.dataSourceA = dataSourceA;
        this.dataSourceB = dataSourceB;
        this.temperature = temperature;
        
        dda = new DataDoubleArray(2);
    }
    
    public DataDoubleArray getData(){
        double[] eAeB = dda.getData();
        
        eAeB[1] = Math.exp(-dataSourceB.getDataAsScalar()/temperature)  
                / Math.exp(-dataSourceA.getDataAsScalar()/temperature);
        eAeB[0] = 1.0;
        
//        System.out.println("meter overlap getdata");
        return dda;
    }
    
    public IDataInfo getDataInfo(){
        return dataInfo;
    }
    public DataTag getTag() {
        return tag;
    }
}
