/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.integrator.IntegratorBox;
import etomica.units.dimensions.Null;

/**
 * Meter used for overlap sampling in the target-sampled system.  The meter
 * measures the ratio of the Boltzmann factors for the reference and target
 * potentials.
 * 
 * @author Tai Boon Tan
 */
public class MeterBoltzmann implements IEtomicaDataSource {
    
    public MeterBoltzmann(IntegratorBox integrator, MeterPotentialEnergy meterPotentialEnergy, double ulatSampled, double ulatMeasured) {
        meterEnergy = new MeterPotentialEnergyFromIntegrator(integrator);
        this.integrator = integrator;
        this.meterPotentialEnergy = meterPotentialEnergy;
        this.ulatSampled = ulatSampled;
        this.ulatMeasured = ulatMeasured;
        data = new DataDoubleArray(2);
        dataInfo = new DataInfoDoubleArray("Scaled Energies", Null.DIMENSION, new int[]{2});
        data.getData()[0] = 1.0;
        tag = new DataTag();
    }

    public IData getData() {
    	double uSampled = meterEnergy.getDataAsScalar() - ulatSampled;
    	double uMeasured = meterPotentialEnergy.getDataAsScalar() - ulatMeasured;
        data.getData()[1] = Math.exp(-(uMeasured - uSampled) / integrator.getTemperature());
        return data;
    }
    
    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    protected final MeterPotentialEnergyFromIntegrator meterEnergy;
    protected final MeterPotentialEnergy meterPotentialEnergy;
    protected final IntegratorBox integrator;
    protected final DataDoubleArray data;
    protected final DataInfoDoubleArray dataInfo;
    protected final double ulatSampled, ulatMeasured;
    protected final DataTag tag;
}
