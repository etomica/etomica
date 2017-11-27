/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.data.DataProcessor;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataInfoFactory;
import etomica.math.function.Function;
import etomica.units.dimensions.Null;

/**
 * DataProcessor that returns the Boltzmann factor of the incoming energy.
 * @author Andrew Schultz
 */
public class BoltzmannProcessor extends DataProcessor {

    private IData data;
    private double temperature;
    private double energyBase;
    
    public IDataInfo processDataInfo(IDataInfo incomingDataInfo) {
        data = incomingDataInfo.makeData();
        IDataInfoFactory factory = incomingDataInfo.getFactory();
        // we get energy in, spit out unitless
        factory.setDimension(Null.DIMENSION);
        dataInfo = factory.makeDataInfo();
        return dataInfo;
    }
    
    public void setEnergyBase(double newEnergyBase) {
        energyBase = newEnergyBase;
    }
    
    public IData processData(IData incomingData) {
        data.E(incomingData);
        data.PE(-energyBase);
        data.TE(-1/temperature);
        data.map(Function.Exp.INSTANCE);
        return data;
    }

    public double getTemperature() {
        return temperature;
    }

    public void setTemperature(double newTemperature) {
        temperature = newTemperature;
    }
}
