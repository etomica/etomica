/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.oneDHardRods;

import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.units.dimensions.Dimension;


public class MeterOverlapSameGaussian1D implements IDataSource {
    DataTag tag;
    DataInfo dataInfo;
    DataSourceScalar dataSourceA, dataSourceB;
    double temperature;
    DataDoubleArray dda;
    
    MeterDifferentImageAdd meterAdd;
    
    
    public static double total, count;
    
    /**
     * Put the system you are measuring in as the first DataSourceScalar
     * @param label
     * @param dimension
     * @param dataSourceA - denominator
     * @param dataSourceB - numerator
     * @param temperature
     */
    public MeterOverlapSameGaussian1D(String label, Dimension dimension, DataSourceScalar dataSourceA,
            DataSourceScalar dataSourceB, double temperature){
        dataInfo = new DataInfoDoubleArray(label, dimension, new int[]{2});
        tag = new DataTag();
        dataInfo.addTag(tag);
        
        
        this.dataSourceA = dataSourceA;
        this.dataSourceB = dataSourceB;
        this.temperature = temperature;
        
        dda = new DataDoubleArray(2);
        
        
        total = 0.0;
        count = 0.0;
    }
    
    public DataDoubleArray getData(){
        double[] eAeB = dda.getData();
        
        double numerator = Math.exp(-dataSourceB.getDataAsScalar()/temperature);
        double gausses = ((MeterDifferentImageAdd1D)dataSourceB).getGaussian();
        double harmonic = 0.5 * (gausses * gausses);
        double denominator = Math.exp(-(dataSourceA.getDataAsScalar()+harmonic)/temperature);
        
        eAeB[1] = numerator / denominator;
        eAeB[0] = 1.0;
     
        
        total += eAeB[1];
        count++;
        
        return dda;
    }
    
    public IDataInfo getDataInfo(){
        return dataInfo;
    }
    public DataTag getTag() {
        return tag;
    }
    
}
