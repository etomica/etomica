/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

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
        mcMoveDiagram = new MCMoveClusterDiagram(integrator.getPotentialMaster());
        integrator.getMoveManager().addMCMove(mcMoveDiagram);
    }

    private static final long serialVersionUID = 1L;
    public MCMoveClusterDiagram mcMoveDiagram;
}
