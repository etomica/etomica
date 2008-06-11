package etomica.modules.chainequilibrium;

import etomica.api.IAtomLeaf;
import etomica.atom.AtomLeafAgentManager;
import etomica.data.DataSourceScalar;
import etomica.units.Fraction;

/**
 * Meter for the reaction conversion (fraction of reacted sites).
 * 
 * @author Andrew Schultz
 */
public class MeterConversion extends DataSourceScalar {

    public MeterConversion(AtomLeafAgentManager agentManager) {
        super("Conversion", Fraction.DIMENSION);
        this.agentManager = agentManager;
        agentIterator = agentManager.makeIterator();
    }

    public double getDataAsScalar() {
        agentIterator.reset();
        int nReacted = 0;
        int total = 0;
        while (agentIterator.hasNext()) {
            IAtomLeaf[] bonds = (IAtomLeaf[])agentIterator.next();
            total += bonds.length;
            for (int i=0; i<bonds.length; i++) {
                if (bonds[i] != null) {
                    nReacted++;
                }
            }
        }
        return ((double)(nReacted))/total;
    }

    protected final AtomLeafAgentManager agentManager;
    protected final AtomLeafAgentManager.AgentIterator agentIterator;
}
