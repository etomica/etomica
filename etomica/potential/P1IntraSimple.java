package etomica.potential;

import etomica.api.IPotentialMaster;
import etomica.atom.iterator.ApiBuilder;


/**
 * Generic intramolecular potential group, having one potential for bonded
 * atoms, and a different potential for unbonded ones.
 *
 * @author David Kofke
 */
public class P1IntraSimple {
    
    public static PotentialGroup makeP1IntraSimple(IPotentialMaster potentialMaster, Potential2 bonded, Potential2 nonbonded) {
        PotentialGroup pGroup = potentialMaster.makePotentialGroup(1);
        pGroup.addPotential(bonded, ApiBuilder.makeAdjacentPairIterator());
        pGroup.addPotential(nonbonded, ApiBuilder.makeNonAdjacentPairIterator());
        return pGroup;
    }
}
   
