/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.oneDHardRods;

import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.units.dimensions.Dimension;


public class MeterOverlap implements IDataSource {
    DataTag tag;
    DataInfo dataInfo;
    DataSourceScalar dataSourceA, dataSourceB;
    double dsABase, dsBBase;
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
    public MeterOverlap(String label, Dimension dimension, DataSourceScalar dataSourceA,
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
        double denominator = -(dataSourceA.getDataAsScalar() - dsABase);
        double power = (numerator - denominator) / temperature;
        eAeB[1] = Math.exp(power);
        
        return dda;
    }
    
    public IDataInfo getDataInfo(){
        return dataInfo;
    }
    public DataTag getTag() {
        return tag;
    }
    public void setDsABase(double dsABase) {
        this.dsABase = dsABase;
    }
    public void setDsBBase(double dsBBase) {
        this.dsBBase = dsBBase;
    }
}
