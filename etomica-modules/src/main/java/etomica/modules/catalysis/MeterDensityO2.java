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
import etomica.units.CompoundDimension;
import etomica.units.Dimension;
import etomica.units.Quantity;
import etomica.units.Volume;

public class MeterDensityO2 extends DataSourceScalar {

    public MeterDensityO2(Box box, ISpecies speciesO, AtomLeafAgentManager interactionAgentManager) {
        super("Density", new CompoundDimension(new Dimension[]{Quantity.DIMENSION, Volume.DIMENSION}, new double[]{+1,-1}));
        this.box = box;
        this.speciesO = speciesO;
        this.interactionAgentManager= interactionAgentManager;
    }

    public double getDataAsScalar() {
        IMoleculeList listO = box.getMoleculeList(speciesO);
        int count = 0;
        for (int i=0; i<listO.getMoleculeCount(); i++) {
            IAtom atom = listO.getMolecule(i).getChildList().getAtom(0);
            CatalysisAgent agent = (CatalysisAgent)interactionAgentManager.getAgent(atom);
            if (agent.isRadical || agent.bondedAtom1.getType() != atom.getType()) continue;
            count++;
        }
        return 0.5*count/box.getBoundary().volume();
    }
    
    private static final long serialVersionUID = 1L;
    protected final Box box;
    protected final ISpecies speciesO;
    protected final AtomLeafAgentManager interactionAgentManager;
}
