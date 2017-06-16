/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.oneDHardRods;

import etomica.data.DataInfo;
import etomica.data.DataSourceScalar;
import etomica.data.DataTag;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.units.Dimension;


public class MeterOverlapSameGaussian implements IEtomicaDataSource {
    DataTag tag;
    DataInfo dataInfo;
    DataSourceScalar dataSourceA, dataSourceB;
    double dsABase, dsBBase;
    double temperature;
    DataDoubleArray dda;
    MeterDifferentImageAdd meterAdd;
    
    /**
     * Put the system you are measuring in as the first DataSourceScalar
     * @param label
     * @param dimension
     * @param dataSourceA - denominator
     * @param dataSourceB - numerator
     * @param temperature
     */
    public MeterOverlapSameGaussian(String label, Dimension dimension, DataSourceScalar dataSourceA,
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
        eAeB[0] = 1.0;
        
        double numerator = -(dataSourceB.getDataAsScalar() - dsBBase);
        //This MUST be done after the call to dataSourceB.getDataAsScalar, in 
        // order to get the right gaussians
        double[] gausses = ((MeterDifferentImageAdd)dataSourceB).getGaussian();
        double harmonic = 0.0;
        for (int i = 0; i < gausses.length; i++){
            harmonic += 0.5 * (gausses[i] * gausses[i]);
        }
        double denominator = -(dataSourceA.getDataAsScalar() + harmonic - dsABase); 
        double power = (numerator - denominator) / temperature;
        eAeB[1] = Math.exp(power); 
        
        return dda;
    }
    
    public IEtomicaDataInfo getDataInfo(){
        return dataInfo;
    }
    public DataTag getTag() {
        return tag;
    }
    public void setDataSourceA(DataSourceScalar dataSourceA) {
        this.dataSourceA = dataSourceA;
    }
    public void setDataSourceB(DataSourceScalar dataSourceB) {
        this.dataSourceB = dataSourceB;
    }

    public void setDsABase(double dsABase) {
        this.dsABase = dsABase;
    }

    public void setDsBBase(double dsBBase) {
        this.dsBBase = dsBBase;
    }
    
}
