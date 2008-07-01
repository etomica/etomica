package etomica.modules.chainequilibrium;

import etomica.api.IAtomLeaf;
import etomica.api.IAtomSet;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.ISpecies;
import etomica.atom.AtomLeafAgentManager;
import etomica.data.DataSourceScalar;
import etomica.units.Fraction;

/**
 * Meter for the reaction conversion (fraction of reacted sites).
 * 
 * @author Andrew Schultz
 */
public class MeterConversion extends DataSourceScalar {

    public MeterConversion(IBox box, AtomLeafAgentManager agentManager) {
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
    
    protected void calcConversion(IAtomSet monomerList) {
        for (int i=0; i<monomerList.getAtomCount(); i++) {
            IAtomLeaf[] bonds = (IAtomLeaf[])agentManager.getAgent(((IMolecule)monomerList.getAtom(i)).getChildList().getAtom(0));
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
    protected final IBox box;
    protected int total, nReacted;
}
