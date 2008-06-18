package etomica.modules.chainequilibrium;

import etomica.api.IAtomSet;
import etomica.api.IAtomTypeLeaf;
import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.api.ISimulation;
import etomica.atom.IAtomKinetic;
import etomica.integrator.IntegratorHard;
import etomica.space.ISpace;

/**
 * IntegratorHard subclass whose sole purpose in life is to allow multiple atoms
 * to be hit with the thermostat at once when using the ANDERSEN_SINGLE
 * thermostat.  By default only one atom is thermostated.
 * 
 * @author Andrew Schultz
 */
public class IntegratorHardThermoFrac extends IntegratorHard {

    public IntegratorHardThermoFrac(ISimulation sim,
            IPotentialMaster potentialMaster, ISpace _space) {
        this(sim, potentialMaster, sim.getRandom(), 0.05, 1.0, _space);
    }

    public IntegratorHardThermoFrac(ISimulation sim,
            IPotentialMaster potentialMaster, IRandom random, double timeStep,
            double temperature, ISpace _space) {
        super(sim, potentialMaster, random, timeStep, temperature, _space);
        setThermostatFrac(0.00001);
    }

    /**
     * Set the fraction of atoms in the box that are thermostated when the
     * ANDERSEN_SINGLE thermostat is used.
     */
    public void setThermostatFrac(double newFrac) {
        thermostatFrac = newFrac;
    }

    /**
     * Returns the fraction of atoms in the box that are thermostated when the
     * ANDERSEN_SINGLE thermostat is used.
     */
    public double getThermostatFrac() {
        return thermostatFrac;
    }

    public void doThermostat() {
        if (--thermostatCount == 0) {
            thermostatCount = thermostatInterval;
            if (thermostat == ThermostatType.ANDERSEN || !initialized) {
                // if initializing the system always randomize the velocity
                randomizeMomenta();
                currentKineticEnergy = meterKE.getDataAsScalar();
            }
            if (thermostat == ThermostatType.VELOCITY_SCALING || !isothermal) {
                scaleMomenta();
                currentKineticEnergy = meterKE.getDataAsScalar();
            }
            else if (thermostat == ThermostatType.ANDERSEN_SINGLE) {
                if (initialized) {
                    IAtomSet atomList = box.getLeafList();
                    int atomCount = atomList.getAtomCount();
                    if (atomCount > 0) {
                        int numThermostatAtoms = (int)Math.ceil(atomCount * thermostatFrac);
                        for (int i=0; i<numThermostatAtoms; i++) {
                            int index = random.nextInt(atomList.getAtomCount());
                            IAtomKinetic a = (IAtomKinetic)atomList.getAtom(index);
                            double m = ((IAtomTypeLeaf)a.getType()).getMass();
                            currentKineticEnergy -= 0.5*m*a.getVelocity().squared();
                            randomizeMomentum(a);
                            currentKineticEnergy += 0.5*m*a.getVelocity().squared();
                        }
                    }
                }
            }
            // ANDERSEN was handled at the start
            else if (thermostat != ThermostatType.ANDERSEN) {
                throw new RuntimeException("Unknown thermostat: "+thermostat);
            }
        }
    }

    protected double thermostatFrac;
}
