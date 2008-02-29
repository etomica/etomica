package etomica.data.meter;

import etomica.EtomicaInfo;
import etomica.data.DataSourceScalar;
import etomica.api.IBox;
import etomica.box.Box;
import etomica.potential.PotentialMaster;
import etomica.units.Energy;

/**
 * Meter for measurement of the total (potential and kinetic) energy in a box
 * This meter is constructed from kinetic-energy and a potential-energy meters
 * An instance of this meter is placed in each box to allow for energy measurements in the box
 */
public final class MeterEnergy extends DataSourceScalar {

    private static final long serialVersionUID = 1L;
    private MeterKineticEnergy kinetic;
    private MeterPotentialEnergy potential;
    
    public MeterEnergy(PotentialMaster potentialMaster) {
    	super("Energy",Energy.DIMENSION);
        kinetic = new MeterKineticEnergy();
        potential = new MeterPotentialEnergy(potentialMaster);
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Total energy (K + P) in a box");
        return info;
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
     * Sets the box(s) where the energy is measured
     * Propagates change to the kinetic- and potential-energy meters
     */
    public void setBox(Box p) {
    	kinetic.setBox(p);
    	potential.setBox(p);
    }
    
    public IBox getBox() {
        return kinetic.getBox();
    }
    
    /**
     * Current value of the total energy (kinetic + potential)
     */
    public double getDataAsScalar() {
        return kinetic.getDataAsScalar() + potential.getDataAsScalar();
    }
}