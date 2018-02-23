/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.catalysis;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.modules.catalysis.InteractionTracker.CatalysisAgent;
import etomica.molecule.IMoleculeList;
import etomica.species.ISpecies;
import etomica.units.dimensions.CompoundDimension;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Quantity;
import etomica.units.dimensions.Volume;

public class MeterDensityCO extends DataSourceScalar {

    public MeterDensityCO(Box box, ISpecies speciesC, AtomLeafAgentManager interactionAgentManager) {
        super("Density", new CompoundDimension(new Dimension[]{Quantity.DIMENSION, Volume.DIMENSION}, new double[]{+1,-1}));
        this.box = box;
        this.speciesC = speciesC;
        this.interactionAgentManager= interactionAgentManager;
    }

    public double getDataAsScalar() {
        IMoleculeList listC = box.getMoleculeList(speciesC);
        int count = 0;
        for (int i = 0; i<listC.size(); i++) {
            IAtom atom = listC.get(i).getChildList().get(0);
            CatalysisAgent agent = (CatalysisAgent)interactionAgentManager.getAgent(atom);
            if (agent.isRadical || agent.bondedAtom2 != null) continue;
            count++;
        }
        return count/box.getBoundary().volume();
    }
    
    private static final long serialVersionUID = 1L;
    protected final Box box;
    protected final ISpecies speciesC;
    protected final AtomLeafAgentManager interactionAgentManager;
}
