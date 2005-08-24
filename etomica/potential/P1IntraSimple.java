package etomica.potential;

import etomica.EtomicaInfo;
import etomica.atom.iterator.ApiBuilder;
import etomica.space.Space;


/**
 * Generic intramolecular potential group, having one potential for bonded
 * atoms, and a different potential for unbonded ones.
 *
 * @author David Kofke
 */
 
 /* History of changes
  * 09/01/02 (DAK) new
  */
 
public class P1IntraSimple extends PotentialGroup implements Potential1.Intramolecular {
    
    public Potential2 bonded;
    public Potential2 nonBonded;
    
    public P1IntraSimple(Space space) {
        super(1, space);
    }
    
    public P1IntraSimple(Space space, Potential2 bonded, Potential2 nonbonded) {
        this(space);
        setBonded(bonded);
        setNonbonded(nonbonded);
    }
    
    /**
     * After constructing bonded potential with this group as its parent, this method
     * must be called with it as an argument to identify it as the bonded potential.
     */
    public void setBonded(Potential2 potential) {
        if(potential == null) return;
        if(bonded != null) removePotential(bonded);
        bonded = potential;
        addPotential(potential, ApiBuilder.makeAdjacentPairIterator());
    }
    
    /**
     * After constructing nonbonded potential with this group as its parent, this method
     * must be called with it as an argument to identify it as the nonbonded potential.
     */
    public void setNonbonded(Potential2 potential) {
        if(potential == null) return;
        if(nonBonded != null) removePotential(nonBonded);
        nonBonded = potential;
        addPotential(potential, ApiBuilder.makeNonAdjacentPairIterator());
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("General intramolecular potential with one bonded and one nonbonded potential");
        return info;
    }

}//end of P1TetheredHardSpheres
   
