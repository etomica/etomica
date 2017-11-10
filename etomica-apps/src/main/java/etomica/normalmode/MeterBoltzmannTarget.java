/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.units.dimensions.Null;

/**
 * Meter used for overlap sampling in the target-sampled system.  The meter
 * measures the ratio of the Boltzmann factors for the harmonic and target
 * potentials.
 * 
 * @author Andrew Schultz
 */
public class MeterBoltzmannTarget implements IDataSource {
    
    public MeterBoltzmannTarget(DataSourceScalar meterTargetEnergy, DataSourceScalar meterRefEnergy) {
        meterEnergy = meterTargetEnergy;
        this.meterRefEnergy = meterRefEnergy;
        data = new DataDoubleArray(2);
        dataInfo = new DataInfoDoubleArray("Scaled Harmonic and hard sphere Energies", Null.DIMENSION, new int[]{2});
        // AccumulatorVirialOverlapSingleAverage expects the first data value to
        // be e_tar / pi_tar, but e_tar = pi_tar
        data.getData()[0] = 1.0;
        tag = new DataTag();
    }

    public void setFrac(double newFrac) {
        frac = newFrac;
    }

    /**
     * Sets the fraction that the system is coupled to measured energy.
     * measured value = exp(-beta*frac*(delta U))
     */
    public IData getData() {
        data.getData()[1] = Math.exp(-frac*(meterRefEnergy.getDataAsScalar() -
                (meterEnergy.getDataAsScalar() - latticeEnergy)) / temperature);
        return data;
    }

    public void setLatticeEnergy(double newLatticeEnergy) {
        latticeEnergy = newLatticeEnergy;
    }
    
    public void setTemperature(double newTemperature) {
        temperature = newTemperature;
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    protected final DataSourceScalar meterEnergy;
    protected final DataSourceScalar meterRefEnergy;
    protected double temperature;
    protected final DataDoubleArray data;
    protected final DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
    protected double latticeEnergy;
    protected double frac;
}
