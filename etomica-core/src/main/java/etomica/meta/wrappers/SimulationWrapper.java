package etomica.meta.wrappers;

import etomica.integrator.Integrator;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorManagerMC;
import etomica.meta.SimulationModel;
import etomica.meta.properties.InstanceProperty;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;

import java.beans.IntrospectionException;
import java.beans.PropertyDescriptor;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Specialized wrapper for Simulations.
 *
 * When a wrapper is created it collects the PotentialMasters in a simulation,
 * which are associated with Integrators internally but conceptually are children of the simulation,
 * and adds a "synthetic" method to return them to the list of the simulation's properties.
 */
public class SimulationWrapper extends ObjectWrapper<Simulation> {

    public SimulationWrapper(Simulation wrapped, SimulationModel simModel, boolean doSerialize) {
        super(wrapped, simModel, doSerialize);

        try {
            PropertyDescriptor descriptor = new PropertyDescriptor(
                    "potentialMasters",
                    SimulationWrapper.class.getDeclaredMethod("getPotentialMasters"),
                    null
            );

            this.properties.add(new InstanceProperty(this, descriptor));
        } catch (IntrospectionException | NoSuchMethodException e) {
            e.printStackTrace();
        }

    }

    public List<PotentialMaster> getPotentialMasters() {
        Set<PotentialMaster> set = new HashSet<>();
        getPotentialMasters(wrapped.getIntegrator(), set);
        return new ArrayList<>(set);
    }

    private static void getPotentialMasters(Integrator integrator, Set<PotentialMaster> set) {
        if(integrator instanceof IntegratorManagerMC) {
            for(Integrator i : ((IntegratorManagerMC) integrator).getIntegrators()) {
                getPotentialMasters(i, set);
            }
        } else if(integrator instanceof IntegratorBox) {
            set.add(((IntegratorBox) integrator).getPotentialMaster());
        }
    }


}
