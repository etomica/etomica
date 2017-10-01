/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.gaussianwork;

import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataSource;
import etomica.data.IDataInfo;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.integrator.IntegratorBox;
import etomica.potential.PotentialMaster;
import etomica.units.dimensions.Null;

/**
 * Multiharmonic system
 * 
 * Meter used for overlap sampling in the A-sampled system. The meter
 * measures the ratio of the Boltzmann factors for System A and System B
 * potentials.
 * 
 * @author Tai Boon Tan
 */
public class MeterBoltzmannA implements IDataSource {
    
    public MeterBoltzmannA(IntegratorBox integratorA, PotentialMaster potentialMasterB) {
    	meterEnergyA = new MeterPotentialEnergyFromIntegrator(integratorA);
    	meterEnergyB = new MeterPotentialEnergy(potentialMasterB);
    	meterEnergyB.setBox(integratorA.getBox());
        data = new DataDoubleArray(2);
        dataInfo = new DataInfoDoubleArray("Scaled System A and System B Energies", Null.DIMENSION, new int[]{2});
    
        data.getData()[0] = 1.0;
        tag = new DataTag();
    }

    public void setTemperature(double newTemperature) {
        temperature = newTemperature;
    }
    
    public double getTemperature() {
        return temperature;
    }
    
    public IData getData() {
        data.getData()[1] = Math.exp(-(meterEnergyB.getDataAsScalar() - 
        		meterEnergyA.getDataAsScalar())/temperature);
        return data;
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    protected MeterPotentialEnergyFromIntegrator meterEnergyA;
    protected MeterPotentialEnergy meterEnergyB;
    protected double temperature;
    protected final DataDoubleArray data;
    protected final DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
}
