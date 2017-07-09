package etomica.meta.wrappers;

import etomica.integrator.Integrator;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorManagerMC;
import etomica.meta.InstanceProperty;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;

import java.beans.IntrospectionException;
import java.beans.PropertyDescriptor;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class SimulationWrapper extends Wrapper<Simulation> {

    public SimulationWrapper(Simulation wrapped) {
        super(wrapped);

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

    public void getPotentialMasters(Integrator integrator, Set<PotentialMaster> set) {
        if(integrator instanceof IntegratorManagerMC) {
            for(Integrator i : ((IntegratorManagerMC) integrator).getIntegrators()) {
                getPotentialMasters(i, set);
            }
        } else if(integrator instanceof IntegratorBox) {
            set.add(((IntegratorBox) integrator).getPotentialMaster());
        }
    }


}
