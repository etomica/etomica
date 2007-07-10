package etomica.potential;

import etomica.atom.AtomAddressManager;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.iterator.AtomIteratorTreeRoot;
import etomica.box.Box;
import etomica.space.Space;

/**
 * Zero-body potential implementing the long-range correction 
 * that compensates for truncation of the potential.  Normally, the
 * concrete instance of this class is defined as an inner class to
 * a PotentialTruncation class.  The PotentialTruncation class defines
 * the scheme used to truncate the potential, and its inner Potential0Lrc
 * subclass implements the formulas needed to integrate the potential
 * and its derivatives over the range affected by the truncation.<br>
 *
 * The components and procedures related to the long-range correction 
 * are as follows.  The Potential2 constructor defines a PotentialTruncation
 * that determines whether and how the pair potential may be truncated.
 * If the PotentialTruncation is not passed as a parameter to the constructor,
 * a default is used as indicated by the Default.TRUNCATE_POTENTIALS boolean
 * flag.  By default, all pair potentials are truncated using the 
 * PotentialTruncationSimple scheme; if TRUNCATE_POTENTIALS is set to
 * false, the PotentialTruncation.Null truncation is applied to all new pair potentials.  The
 * PotentialTruncation for a potential cannot be changed after the potential
 * is instantiated.<br>
 *
 * Each PotentialTruncation class defines an inner Potential0Lrc subclass that 
 * provides the appropriate methods for computing the long-range correction 
 * to the energy and its first two derivatives.  This class is instantiated
 * by the PotentialTruncation class in its constructor.  Upon its instantiation,
 * the Potential0Lrc class is added to the group of long-range correction potentials
 * that is kept by a single Potential0GroupLrc instance in the PotentialMaster.<br>
 *
 * Before the calculate method of PotentialMaster is called to compute something,
 * its set(Box) method must have been previously called, which identifies to
 * all potentials (including Potential0GroupLrc) which box is subject to the 
 * ensuing calculation.  Potential0Group ignores this notification if the
 * given box is the same as the one specified in the previous call; otherwise
 * it passes the identified box to all the set(Box) methods (inherited from Potential0)
 * of the Potential0Lrc classes it holds.<br>
 *
 * Then when the calculate(IteratorDirective, PotentialCalculation) method of
 * PotentialMaster is invoked, it passes the call on to the calculate methods of
 * its child potentials, including the Potential0GroupLrc instance if present.
 * Potential0GroupLrc checks that a box has been specified, that its
 * enableLrc flag is <b>true</b> (the default), and that the given iteratorDirective's
 * includeP0Lrc flag is also <b>true</b> (default is <b>false</b>).  If so, it 
 * calls the calculate methods of all child Potential0Lrc classes.
 *
 * The Potential0Lrc class will use the volume from the specified box and the
 * size method of the iterator of its associated potential to determine the
 * pair density.  The Potential0Lrc methods are called if the PotentialCalculation
 * implements Potential0Calculation.
 *
 * @author David Kofke
 */
public abstract class Potential0Lrc extends Potential0 implements PotentialSoft {
    
    protected final AtomType[] types;
    protected final boolean interType;
    protected final Potential truncatedPotential;
    protected final int[] lrcAtomsPerMolecule = new int[2];
//    protected final ISpeciesAgent[] agents = new ISpeciesAgent[2];
    
    public Potential0Lrc(Space space, AtomType[] types, Potential truncatedPotential) {
        super(space);
        this.types = (AtomType[])types.clone();
        if(types.length != 2) {
            throw new IllegalArgumentException("LRC developed only for two-body potentials; must give two species to constructor");
        }
        interType = (types[0] != types[1]);
        this.truncatedPotential = truncatedPotential;
        divisor = 1;
    }  
    
    /**
     * Returns the potential whose truncation this lrcPotential exists to correct. 
     */
    public Potential getTruncatedPotential() {
        return truncatedPotential;
    }

    public void setBox(Box p) {
        if (lrcAtomsPerMolecule[0] == 0 || lrcAtomsPerMolecule[1] == 0) {
            // count the number of Atoms of the relevant type in each molecule
            AtomIteratorTreeRoot treeIterator = new AtomIteratorTreeRoot();
            for (int i=0; i<2; i++) {
                if (lrcAtomsPerMolecule[i] == 0) {
                    if (types[i].getDepth() == AtomAddressManager.MOLECULE_DEPTH) {
                        lrcAtomsPerMolecule[i] = 1;
                    }
                    else if (p.getNMolecules(types[i].getSpecies()) > 0) {
                        treeIterator.setRootAtom(p.getMoleculeList(types[i].getSpecies()).getAtom(0));
                        treeIterator.setIterationDepth(types[i].getDepth());
                        treeIterator.reset();
                        int numAtoms = 0;
                        for (IAtom atom = treeIterator.nextAtom(); atom != null;
                             atom = treeIterator.nextAtom()) {
                            if (atom.getType() == types[i]) numAtoms++;
                        }
                        lrcAtomsPerMolecule[i] = numAtoms;
                    }
                }
            }
        }
        box = p;
    }
    
    public void setTargetAtoms(IAtom targetAtom) {
        if (targetAtom == null) {
            divisor = 1;
            return;
        }
        int typeIndex = 1;
        if (types[0].isDescendedFrom(targetAtom.getType())) {
            typeIndex = 0;
        }
        else if (!types[1].isDescendedFrom(targetAtom.getType())) {
            divisor = 0;
            return;
        }
        divisor = box.getNMolecules(types[typeIndex].getSpecies()) * lrcAtomsPerMolecule[typeIndex];
        if (!interType) {
            divisor = (divisor - 1)/2.0;
        }
    }
     
    /**
     * Returns the number of pairs formed from molecules of the current
     * species, in the given box.
     * @param box
     * @return int
     */
    protected int nPairs() {
        if(divisor == 0) return 0;
        int nPairs = 0;
        int n0 = box.getNMolecules(types[0].getSpecies())*lrcAtomsPerMolecule[0];
        if(interType) {
            int n1 = box.getNMolecules(types[1].getSpecies());
            nPairs = n0*n1*lrcAtomsPerMolecule[1];
        }
        else {
            nPairs = n0*(n0-1)/2;
        }
        return (int)Math.round(nPairs/divisor);
    }
    
    /**
     * Long-range correction to the energy u.
     */
    public abstract double uCorrection(double pairDensity);
        
    /**
    * Long-range correction to r*du/dr.
    */
    public abstract double duCorrection(double pairDensity);
    
    /**
    * Long-range correction to r^2 d2u/dr2.
    */
    public abstract double d2uCorrection(double pairDensity);
    
    protected double divisor;
    protected Box box;

}//end of Potential0Lrc
