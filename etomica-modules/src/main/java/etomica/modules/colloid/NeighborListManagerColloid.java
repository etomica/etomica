/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.colloid;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.molecule.IMoleculeList;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.IPotentialAtomic;
import etomica.species.ISpecies;

public class NeighborListManagerColloid extends NeighborListManager {

    public NeighborListManagerColloid(PotentialMasterList potentialMasterList,
                                      double range, Box box) {
        super(potentialMasterList, range, box);
    }
    
    public void setSpeciesColloid(ISpecies newSpeciesColloid) {
        speciesColloid = newSpeciesColloid;
    }
    
    public void setPotentialMC(IPotentialAtomic newPotentailMC) {
        p2mc = newPotentailMC;
    }
    
    public void setSpeciesMonomer(ISpecies newSpeciesMonomer) {
        speciesMonomer = newSpeciesMonomer;
    }
    
    public void setPotentialPseudo(IPotentialAtomic newPotentailPseudo) {
        p2pseudo = newPotentailPseudo;
    }
    
    public void setChainLength(int newChainLength) {
        chainLength = newChainLength;
    }
    
    protected void neighborSetup() {
        super.neighborSetup();
        
        IAtom colloidAtom = box.getMoleculeList(speciesColloid).get(0).getChildList().get(0);
        IMoleculeList monomers = box.getMoleculeList(speciesMonomer);
        if (monomers.size() == 0) {
            return;
        }

        int colloidTypeIndex = colloidAtom.getType().getIndex();
        int monomerTypeIndex = monomers.get(0).getChildList().get(0).getType().getIndex();
        for (int i = 0; i<monomers.size(); i++) {
            IAtom atom = monomers.get(i).getChildList().get(0);
            agentManager2Body.getAgent(colloidAtom).addUpNbr(atom, colloidTypeIndex);
            agentManager2Body.getAgent(atom).addDownNbr(colloidAtom, monomerTypeIndex);
        }
        
        for (int i = 0; i<monomers.size(); i+=chainLength) {
            IAtom iAtom = monomers.get(i).getChildList().get(0);
            for (int j = i+chainLength; j<monomers.size(); j+=chainLength) {
                IAtom jAtom = monomers.get(j).getChildList().get(0);
                // iAtom and jAtom are both bonded to the colloid.  we need to keep them apart with p2seudo
                agentManager2Body.getAgent(iAtom).addUpNbr(jAtom, monomerTypeIndex);
                agentManager2Body.getAgent(jAtom).addDownNbr(iAtom, monomerTypeIndex);
            }
        }
    }

    protected ISpecies speciesColloid, speciesMonomer;
    protected IPotentialAtomic p2mc, p2pseudo;
    protected int chainLength;

    public static class NeighborListAgentSourceColloid extends PotentialMasterList.NeighborListAgentSource {
        public NeighborListAgentSourceColloid(double range) {
            super(range);
        }

        public NeighborListManagerColloid makeAgent(Box box) {
            return new NeighborListManagerColloid(potentialMaster, range, box);
        }
    }
}
