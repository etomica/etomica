/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.box.Box;
import etomica.potential.PotentialMaster;
import etomica.data.DataSourceScalar;
import etomica.units.dimensions.Energy;

/**
 * Meter for measurement of the total (potential and kinetic) energy in a box
 * This meter is constructed from kinetic-energy and a potential-energy meters
 * An instance of this meter is placed in each box to allow for energy measurements in the box
 */
public final class MeterEnergy extends DataSourceScalar {

    private static final long serialVersionUID = 1L;
    protected DataSourceScalar kinetic;
    protected DataSourceScalar potential;
    
    public MeterEnergy(PotentialMaster potentialMaster, Box box) {
    	super("Energy",Energy.DIMENSION);
        kinetic = new MeterKineticEnergy();
        ((MeterKineticEnergy)kinetic).setBox(box);
        potential = new MeterPotentialEnergy(potentialMaster);
        ((MeterPotentialEnergy)potential).setBox(box);
    }
    
    public DataSourceScalar getKinetic() {
        return kinetic;
    }

    public void setKinetic(DataSourceScalar kinetic) {
        this.kinetic = kinetic;
    }

    public DataSourceScalar getPotential() {
        return potential;
    }

    public void setPotential(DataSourceScalar potential) {
        this.potential = potential;
    }

    /**
     * @return the current value of the total kinetic energy of the molecules in the box
     */
    public double getKineticEnergy() {return kinetic.getDataAsScalar();}
    /**
     * @return the current value of the total potential energy of the molecules in the box
     */
    public double getPotentialEnergy() {return potential.getDataAsScalar();}
     
    /**
     * Current value of the total energy (kinetic + potential)
     */
    public double getDataAsScalar() {
        return kinetic.getDataAsScalar() + potential.getDataAsScalar();
    }
}
