package etomica;

import etomica.units.Unit;
import etomica.units.Dimension;

/**
 * Meter for measurement of the total (potential and kinetic) energy in a phase
 * This meter is constructed from kinetic-energy and a potential-energy meters
 * An instance of this meter is placed in each phase to allow for energy measurements in the phase
 */
public final class MeterEnergy extends Meter implements EtomicaElement
{
    private MeterKineticEnergy kinetic;
    private MeterPotentialEnergy potential;
    
    public MeterEnergy() {
        this(Simulation.instance);
    }
    public MeterEnergy(Simulation sim) 
    {
        super(sim);
        setLabel("Energy");
        kinetic = new MeterKineticEnergy();
        potential = new MeterPotentialEnergy();
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Total energy (K + P) in a phase");
        return info;
    }

    /**
     * Declaration whether this meter uses the boundary object of phase when making its measurements
     */
    public final boolean usesPhaseBoundary() {
        return kinetic.usesPhaseBoundary() || potential.usesPhaseBoundary();
    }
    /**
     * Declaration whether this meter uses the iteratorFactory of phase when making its measurements
     */
    public final boolean usesPhaseIteratorFactory() {
        return kinetic.usesPhaseIteratorFactory() || potential.usesPhaseIteratorFactory();
    }

    public Dimension getDimension() {return Dimension.ENERGY;}

    /**
     * @return a handle to the kinetic-energy meter
     */
    public MeterKineticEnergy meterKinetic() {return kinetic;}
    /**
     * @return a handle to the potential-energy meter
     */
    public MeterPotentialEnergy meterPotential() {return potential;}
    /**
     * Accessor method to set the kinetic-energy meter to something other than the default
     */
    public void setMeterKinetic(MeterKineticEnergy mke) {kinetic = mke;}
    /**
     * Accessor method to set the potential-energy meter to something other than the default
     */
    public void setMeterPotential(MeterPotentialEnergy mpe) {potential = mpe;}
    
    /**
     * @return the current value of the total kinetic energy of the molecules in the phase
     */
    public double kinetic() {return kinetic.currentValue();}
    /**
     * @return the current value of the total potential energy of the molecules in the phase
     */
    public double potential() {return potential.currentValue();}
    
    /**
     * Sets the phase where the energy is measured
     * Propagates change to the kinetic- and potential-energy meters
     */
    public void setPhase(Phase p) {
        super.setPhase(p);
        kinetic.setPhase(p);
        potential.setPhase(p);
    }
    
    /**
     * Current value of the total energy (kinetic + potential)
     */
    public double currentValue()
    {
        return kinetic.currentValue() + potential.currentValue();
    }
}