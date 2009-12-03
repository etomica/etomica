package etomica.models.oneDHardRods;

import etomica.api.IRandom;
import etomica.data.DataInfo;
import etomica.data.DataSourceScalar;
import etomica.data.DataTag;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.units.Dimension;


public class MeterOverlapSameGaussian implements IEtomicaDataSource {
    DataTag tag;
    DataInfo dataInfo;
    DataSourceScalar dataSourceA, dataSourceB;
    double temperature;
    DataDoubleArray dda;
    
    MeterDifferentImageAdd meterAdd;
    double[] gaussian;
    protected final IRandom random;
    
    /**
     * Put the system you are measuring in as the first DataSourceScalar
     * @param label
     * @param dimension
     * @param dataSourceA - denominator
     * @param dataSourceB - numerator
     * @param temperature
     */
    MeterOverlapSameGaussian(String label, Dimension dimension, DataSourceScalar dataSourceA,
            DataSourceScalar dataSourceB, double temperature, IRandom rand){
        dataInfo = new DataInfoDoubleArray(label, dimension, new int[]{2});
        tag = new DataTag();
        dataInfo.addTag(tag);
        
        this.random = rand;
        this.dataSourceA = dataSourceA;
        this.dataSourceB = dataSourceB;
        this.temperature = temperature;
        
        dda = new DataDoubleArray(2);
        
        gaussian = new double[2];
    }
    
    public DataDoubleArray getData(){
        double[] eAeB = dda.getData();
        
        generateGaussian();
        
        
        
        double numerator = Math.exp(-dataSourceB.getDataAsScalar()/temperature);
        double denominator = Math.exp(-dataSourceA.getDataAsScalar()/temperature);
        
        eAeB[1] = numerator / denominator;
        eAeB[0] = 1.0;
        
        return dda;
    }
    
    public IEtomicaDataInfo getDataInfo(){
        return dataInfo;
    }
    public DataTag getTag() {
        return tag;
    }
    
    public void generateGaussian(){
        double sqrtT = Math.sqrt(temperature);
        double realGauss = random.nextGaussian() * sqrtT;
        double imagGauss = random.nextGaussian() * sqrtT;

        gaussian[0] = realGauss;
        gaussian[1] = imagGauss;
    }
    
    
    

}
