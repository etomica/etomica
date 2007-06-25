package etomica.virial.simulations;

import etomica.space.Space;
import etomica.virial.ClusterAbstract;
import etomica.virial.ClusterWeight;
import etomica.virial.MCMoveClusterDiagram;
import etomica.virial.SpeciesFactory;

/**
 * Virial simulation using diagram sampling.
 *
 * @author Andrew Schultz
 */
public class SimulationVirialDiagram extends SimulationVirial {

    public SimulationVirialDiagram(Space space,
            SpeciesFactory speciesFactory, double temperature,
            ClusterWeight aSampleCluster, ClusterAbstract refCluster,
            ClusterAbstract[] targetClusters) {
        super(space, speciesFactory, temperature, aSampleCluster,
                refCluster, targetClusters);
        mcMoveDiagram = new MCMoveClusterDiagram(integrator.getPotential());
        integrator.getMoveManager().addMCMove(mcMoveDiagram);
    }

    private static final long serialVersionUID = 1L;
    public MCMoveClusterDiagram mcMoveDiagram;
}
