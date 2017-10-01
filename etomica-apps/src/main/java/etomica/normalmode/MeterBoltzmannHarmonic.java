/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataSource;
import etomica.data.IEtomicaDataInfo;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.potential.PotentialMaster;
import etomica.units.dimensions.Null;

/**
 * Meter used for overlap sampling in the harmonic-sampled system.  The meter
 * measures the ratio of the Boltzmann factors for the harmonic and target
 * potentials.
 * 
 * @author Andrew Schultz
 */
public class MeterBoltzmannHarmonic implements IDataSource {
    
    public MeterBoltzmannHarmonic(MCMoveHarmonic mcMoveHarmonic, PotentialMaster potentialMaster) {
        this.mcMoveHarmonic = mcMoveHarmonic;
        meterEnergy = new MeterPotentialEnergy(potentialMaster);
        meterEnergy.setBox(mcMoveHarmonic.getBox());
        data = new DataDoubleArray(2);
        dataInfo = new DataInfoDoubleArray("Scaled Harmonic and hard sphere Energies", Null.DIMENSION, new int[]{2});
        // AccumulatorVirialOverlapSingleAverage expects the first data value to
        // be e_ref / pi_ref, but e_ref = pi_ref
        data.getData()[0] = 1.0;
        tag = new DataTag();
    }

    /**
     * Sets the fraction that the system is coupled to measured energy.
     * measured value = exp(-beta*frac*(delta U))
     */
    public void setFrac(double newFrac) {
        frac = newFrac;
    }

    public void setTemperature(double newTemperature) {
        temperature = newTemperature;
    }
    
    public double getTemperature() {
        return temperature;
    }
    
    public IData getData() {
        data.getData()[1] = Math.exp(-frac*((meterEnergy.getDataAsScalar() - latticeEnergy) - 
                mcMoveHarmonic.getLastTotalEnergy())/temperature);
        return data;
    }

    public void setLatticeEnergy(double newLatticeEnergy) {
        latticeEnergy = newLatticeEnergy;
    }
    
    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    protected double temperature;
    protected final MeterPotentialEnergy meterEnergy;
    protected final MCMoveHarmonic mcMoveHarmonic;
    protected final DataDoubleArray data;
    protected final DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
    protected double latticeEnergy;
    protected double frac;
}
