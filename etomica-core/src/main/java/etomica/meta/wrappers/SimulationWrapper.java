package etomica.meta.wrappers;

import etomica.atom.IAtomList;
import etomica.integrator.Integrator;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorManagerMC;
import etomica.meta.SimulationModel;
import etomica.meta.properties.ArrayProperty;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;

import java.beans.IntrospectionException;
import java.beans.PropertyDescriptor;
import java.util.HashSet;
import java.util.Set;

/**
 * Specialized wrapper for Simulations.
 *
 * When a wrapper is created it collects the PotentialMasters in a simulation,
 * which are associated with Integrators internally but conceptually are children of the simulation,
 * and adds a "synthetic" method to return them to the list of the simulation's properties.
 */
public class SimulationWrapper extends ObjectWrapper<Simulation> {

    private double boxes[][][];

    public SimulationWrapper(Simulation wrapped, SimulationModel simModel, boolean doSerialize) {
        super(wrapped, simModel, doSerialize);

        try {
            PropertyDescriptor descriptor = new PropertyDescriptor(
                    "potentialMasters",
                    SimulationWrapper.class.getDeclaredMethod("getPotentialMasters"),
                    null
            );

            this.childProps.add(new ArrayProperty(this, descriptor));
        } catch (IntrospectionException | NoSuchMethodException e) {
            e.printStackTrace();
        }
    }

    public PotentialMaster[] getPotentialMasters() {
        Set<PotentialMaster> set = new HashSet<>();
        getPotentialMasters(wrapped.getIntegrator(), set);
        return set.toArray(new PotentialMaster[set.size()]);
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

    /**
     * Array of coordinates of all atoms in all boxes in the wrapped Simulation.
     * @return 3-dimensional array of coordinates. First index indicates the box; second index
     * indicates the leaf atom in the box; third index indicates the coordinate (i.e., x, y, z)
     */
    public double[][][] getAllCoordinates() {
        Simulation sim = simModel.getSimulation();
        if(boxes == null || boxes.length != sim.getBoxCount()) {
            boxes = new double[sim.getBoxCount()][][];
        }

        for(int i = 0; i < sim.getBoxCount(); i++) {
            IAtomList leafList = sim.getBox(i).getLeafList();
            if(boxes[i] == null || boxes[i].length != leafList.size()) {
                boxes[i] = new double[leafList.size()][sim.getSpace().D()];
            }

            for(int j = 0; j < leafList.size(); j++) {
                leafList.get(j).getPosition().assignTo(boxes[i][j]);
            }
        }
        return boxes;
    }


/*
    public static void main(String[] args) {
        Simulation sim = new HSMD2D();
        SimulationWrapper wrapper = new SimulationWrapper(sim, new SimulationModel(sim), true);
        double[][][] coords = wrapper.getAllCoordinates();
        System.out.println("done");

    }
*/
}
