package etomica.virial;

import etomica.atom.AtomHydrogen;
import etomica.atom.AtomTypeOriented;
import etomica.chem.elements.Hydrogen;
import etomica.config.ConformationLinear;
import etomica.integrator.IntegratorMC;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;
import etomica.units.Kelvin;
import etomica.virial.cluster.ClusterAbstract;
import etomica.virial.cluster.ClusterWeight;
import etomica.virial.mcmove.MCMoveClusterRingRegrowOrientation;
import org.junit.jupiter.api.Test;

class MCMoveClusterRingRegrowOrientationTest {

    @Test
    public void testStability() {
        Space space = Space3D.getInstance();
        ClusterWeight cluster = new ClusterWeight() {

            @Override
            public double value(BoxCluster box) {
                return 1;
            }

            @Override
            public void setTemperature(double temperature) {
            }

            @Override
            public int pointCount() {
                return 1;
            }

            @Override
            public ClusterAbstract makeCopy() {
                return null;
            }
        };
        for (int p = 2; p <= 512; p *= 2) {
            Simulation sim = new Simulation(space);
            SpeciesGeneral species = new SpeciesBuilder(space)
                    .addCount(new AtomTypeOriented(Hydrogen.INSTANCE, space.makeVector()), p)
                    .withConformation(new ConformationLinear(space))
                    .withAtomFactory(atomType -> {
                        return new AtomHydrogen(space, ((AtomTypeOriented) atomType), 1);
                    })
                    .build();
            sim.addSpecies(species);
            BoxCluster box = new BoxCluster(cluster, space);
            sim.addBox(box);
            box.setNMolecules(species, 2);
            IntegratorMC integrator = new IntegratorMC(new PotentialComputeAggregate(), sim.getRandom(), 1, box);
            MCMoveClusterRingRegrowOrientation move = new MCMoveClusterRingRegrowOrientation(sim.getRandom(), space, p);

            for (int iTemp = 40; iTemp <= 40; iTemp += 2) {
                move.setStiffness(Kelvin.UNIT.toSim(iTemp), species.getAtomType(0).getMass());
                integrator.getMoveManager().addMCMove(move);
                integrator.reset();
                int total = 100;
                for (int i = 0; i < total; i++) {
                    integrator.doStep();
                }
                System.out.println("p = " + p + " ,Temp = " + iTemp + " ,acceptance ratio = " + move.getTracker().acceptanceRatio());
            }

        }
    }

}