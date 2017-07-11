/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.chainequilibrium;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.molecule.IMoleculeList;
import etomica.species.ISpecies;
import etomica.units.dimensions.Fraction;

/**
 * Meter for the reaction conversion (fraction of reacted sites).
 * 
 * @author Andrew Schultz
 */
public class MeterConversion extends DataSourceScalar {

    public MeterConversion(Box box, AtomLeafAgentManager agentManager) {
        super("Conversion", Fraction.DIMENSION);
        this.box = box;
        this.agentManager = agentManager;
    }

    public void setSpecies(ISpecies[] newSpecies) {
        species = newSpecies;
    }

    public double getDataAsScalar() {
        total = 0;
        nReacted = 0;
        if (species == null) {
            calcConversion(box.getMoleculeList());
        }
        else {
            for (int i=0; i<species.length; i++) {
                calcConversion(box.getMoleculeList(species[i]));
            }
        }
        return ((double)(nReacted))/total;
    }
    
    protected void calcConversion(IMoleculeList monomerList) {
        for (int i=0; i<monomerList.getMoleculeCount(); i++) {
            IAtom[] bonds = (IAtom[])agentManager.getAgent(monomerList.getMolecule(i).getChildList().getAtom(0));
            total += bonds.length;
            for (int j=0; j<bonds.length; j++) {
                if (bonds[j] != null) {
                    nReacted++;
                }
            }
        }
    }

    protected final AtomLeafAgentManager agentManager;
    protected ISpecies[] species;
    protected final Box box;
    protected int total, nReacted;
}
