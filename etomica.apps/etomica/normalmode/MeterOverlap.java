package etomica.normalmode;

import etomica.data.DataTag;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.units.Null;


public class MeterOverlap implements IEtomicaDataSource {
    protected DataTag tag;
    protected IEtomicaDataInfo dataInfo;
    protected IEtomicaDataSource dataSourceSame, dataSourceDifferent;
    protected double temperature;
    protected DataDoubleArray dda;
    
    /**
     * Meter to measure the "overlap" function for free energy calculations,
     * used as input to an AccumulatorVirialOverlapSingleAverage.
     * 
     * @param dataSourceSame - data source that returns the energy of the
     *    system being sampled
     * @param dataSourceDifferent - data source that returns the energy of the
     *    system not being sampled
     * @param temperature - the temperature
     */
    public MeterOverlap(IEtomicaDataSource dataSourceSame,
            IEtomicaDataSource dataSourceDifferent, double temperature){
        dataInfo = new DataInfoDoubleArray("overlap", Null.DIMENSION, new int[]{2});
        tag = new DataTag();
        dataInfo.addTag(tag);
        
        this.dataSourceSame = dataSourceSame;
        this.dataSourceDifferent = dataSourceDifferent;
        this.temperature = temperature;
        
        dda = new DataDoubleArray(2);
        double[] eAeB = dda.getData();
        eAeB[0] = 1.0;
    }
    
    public DataDoubleArray getData(){
        double[] eAeB = dda.getData();

        eAeB[1] = Math.exp((dataSourceSame.getData().getValue(0)-dataSourceDifferent.getData().getValue(0))/temperature);  
        return dda;
    }
    
    public IEtomicaDataInfo getDataInfo(){
        return dataInfo;
    }
    public DataTag getTag() {
        return tag;
    }
}
