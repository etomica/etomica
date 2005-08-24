package etomica.data.meter;

import etomica.EtomicaInfo;
import etomica.Phase;
import etomica.data.DataSourceScalar;
import etomica.potential.PotentialMaster;
import etomica.units.Dimension;

/**
 * Meter for measurement of the total (potential and kinetic) energy in a phase
 * This meter is constructed from kinetic-energy and a potential-energy meters
 * An instance of this meter is placed in each phase to allow for energy measurements in the phase
 */
public final class MeterEnergy extends DataSourceScalar implements Meter {

    private MeterKineticEnergy kinetic;
    private MeterPotentialEnergy potential;
    
    public MeterEnergy(PotentialMaster potentialMaster) {
    	super("Energy",Dimension.ENERGY);
        kinetic = new MeterKineticEnergy();
        potential = new MeterPotentialEnergy(potentialMaster);
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Total energy (K + P) in a phase");
        return info;
    }
    
    /**
     * @return the current value of the total kinetic energy of the molecules in the phase
     */
    public double getKineticEnergy() {return kinetic.getDataAsScalar();}
    /**
     * @return the current value of the total potential energy of the molecules in the phase
     */
    public double getPotentialEnergy() {return potential.getDataAsScalar();}
    
    /**
     * Sets the phase(s) where the energy is measured
     * Propagates change to the kinetic- and potential-energy meters
     */
    public void setPhase(Phase p) {
    	kinetic.setPhase(p);
    	potential.setPhase(p);
    }
    
    public Phase getPhase() {
        return kinetic.getPhase();
    }
    
    /**
     * Current value of the total energy (kinetic + potential)
     */
    public double getDataAsScalar() {
        return kinetic.getDataAsScalar() + potential.getDataAsScalar();
    }
}