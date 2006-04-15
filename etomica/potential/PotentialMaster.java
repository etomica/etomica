package etomica.potential;

import java.util.Arrays;

import etomica.atom.Atom;
import etomica.atom.AtomSet;
import etomica.atom.AtomType;
import etomica.atom.iterator.AtomIteratorAll;
import etomica.atom.iterator.AtomsetIteratorPDT;
import etomica.atom.iterator.AtomsetIteratorSinglet;
import etomica.atom.iterator.IteratorDirective;
import etomica.atom.iterator.IteratorFactory;
import etomica.phase.Phase;
import etomica.space.Space;
import etomica.species.Species;


/**
 * Manager of all potentials in simulation.
 * Most calls to compute the energy or other potential calculations begin
 * with the calculate method of this class.  It then passes the calculation 
 * on to the contained potentials.
 *
 * @author David Kofke
 */
 
public class PotentialMaster implements java.io.Serializable {
    
    public PotentialMaster(Space space) {
        this(space,IteratorFactory.INSTANCE);
    } 
    
    public PotentialMaster(Space space, IteratorFactory iteratorFactory) {
        this.space = space;
        this.iteratorFactory = iteratorFactory;
    }
    
	/**
	 * Returns the object that oversees the long-range
	 * correction zero-body potentials.
	 */
	 public PotentialMasterLrc lrcMaster() {
		if(lrcMaster == null) lrcMaster = new PotentialMasterLrc(space);
		return lrcMaster;
	 }
     
     public PotentialGroup makePotentialGroup(int nBody) {
         return new PotentialGroup(nBody,space);
     }

     /**
      * Performs the given PotentialCalculation on the atoms of the given Phase.
      * Sets the phase for all molecule iterators and potentials, sets target
      * and direction for iterators as specified by given IteratorDirective,
      * and applies doCalculation of given PotentialCalculation with the iterators
      * and potentials.
      */
    public void calculate(Phase phase, IteratorDirective id, PotentialCalculation pc) {
        if(!enabled) return;
    	Atom targetAtom = id.getTargetAtom();
    	boolean phaseChanged = (phase != mostRecentPhase);
    	mostRecentPhase = phase;
    	for(PotentialLinker link=first; link!=null; link=link.next) {
    	    if(!link.enabled) continue;
    	    if(phaseChanged) {
    	        link.iterator.setPhase(phase);
    	        link.potential.setPhase(phase);
    	    }
    	    link.iterator.setTarget(targetAtom);
    	    link.iterator.setDirection(id.direction());
    	    pc.doCalculation(link.iterator, id, link.potential);
        }//end for
        if(lrcMaster != null) {
            lrcMaster.calculate(phase, id, pc);
        }
    }//end calculate
    
    /**
     * Indicates to the PotentialMaster that the given potential should apply to 
     * the specified species.  Exception is thrown if the potential.nBody() value
     * is different from the length of the species array.  Thus, for example, if
     * giving a 2-body potential, then the array should contain exactly
     * two species; the species may refer to the same instance (appropriate for an 
     * intra-species potential, defining the iteractions between molecules of the
     * same species).
     */
    public void addPotential(Potential potential, Species[] species) {
    	if (potential.nBody() == 0) {
    		addPotential(potential, new AtomIterator0(),null);
    	}
        else if (potential.nBody() == Integer.MAX_VALUE) {
            addPotential(potential, new AtomIteratorAll(species), null);
        }
    	else if (potential.nBody() != species.length) {
    		throw new IllegalArgumentException("Illegal species length");
    	}
        else {
            AtomsetIteratorPDT iterator = iteratorFactory.makeMoleculeIterator(species);
            addPotential(potential, iterator, moleculeTypes(species));
            if(potential instanceof PotentialTruncated) {
                Potential0Lrc lrcPotential = ((PotentialTruncated)potential).makeLrcPotential(moleculeTypes(species)); 
                if(lrcPotential != null) {
                    lrcMaster().addPotential(
                        lrcPotential,
                        new AtomIterator0(),null);
                }
            }
        }
    }
    
    /**
     * Indicates to the PotentialMaster that the given potential should apply to 
     * the specified atom types.  The potential is assumed to be intermolecular.
     * The given types should not include any type which is the descendent of 
     * another.  Potential group hierarchy will be constructed as needed above
     * the level of the given atom types.
     * <p>
     * The order of the elements in the atomTypes array is not relevant, and is
     * subject to rearrangement by the method -- the array is sorted (using the compareTo
     * method of AtomType) before doing anything else.
     * 
     */
    public void addPotential(Potential potential, AtomType[] atomTypes) {
        if (potential.nBody() != atomTypes.length) {
            throw new IllegalArgumentException("nBody of potential must match number of atom types");
        }
        Arrays.sort(atomTypes);
        // depth of molecules
        int maxDepth = 3;
        int[] depth = new int[atomTypes.length];
        for (int i=0; i<atomTypes.length; i++) {
            depth[i] = atomTypes[i].getAddressManager().getDepth();
            if (depth[i] > maxDepth) maxDepth = depth[i];
        }
        if (maxDepth == 3) {
            addPotential(potential,moleculeSpecies(atomTypes));
            return;
        }
        AtomType[] parentAtomTypes = new AtomType[atomTypes.length];
        for (int i=0; i<atomTypes.length; i++) {
            if (depth[i] == maxDepth) {
                parentAtomTypes[i] = atomTypes[i].getParentType();
            }
            else {
                parentAtomTypes[i] = atomTypes[i];
            }
        }
        // look for a PotentialGroup that applies to parentAtomTypes
        PotentialGroup pGroup = getPotential(parentAtomTypes);
        if (pGroup == null) { // didn't find an appropriate potentialgroup
            pGroup = makePotentialGroup(atomTypes.length);
            addPotential(pGroup,parentAtomTypes);
        }
        pGroup.addPotential(potential,atomTypes);
    }
    
    /**
     * Notifies the PotentialMaster that the sub-potential has been added to 
     * the given PotentialGroup, which is associated (but not necessarily held 
     * by) this PotentialMaster.
     * This method is called by PotentialGroup and should not be called in
     * other circumstances.
     */
    public void potentialAddedNotify(Potential subPotential, PotentialGroup pGroup) {
        // do nothing.  this is here for subclasses to override
    }
    
    /**
     * Returns the potential that applies to the specified types,
     * or null of no existing potential applies.
     */
    public PotentialGroup getPotential(AtomType[] types) {
        for(PotentialLinker link=first; link!=null; link=link.next) {
            if (link.potential instanceof PotentialGroup) {
                if(Arrays.equals(types,link.types)) {
                    return (PotentialGroup)link.potential;
                }
                PotentialGroup candidate = ((PotentialGroup)link.potential).getPotential(types);
                if (candidate != null) {
                    return candidate;
                }
            }
        }
        return null;
    }
    
    /**
     * Returns the AtomTypes that the given potential applies to if the given 
     * potential is within this potential group.  If the potential is not 
     * contained by the potential master or any PotentialGroup it holds, or 
     * does not apply to specific AtomTypes, null is returned.
     */
    public AtomType[] getAtomTypes(Potential potential) {
        for(PotentialLinker link=first; link!=null; link=link.next) {
            if (link.potential == potential) {
                return link.types;
            }
            if (link.potential instanceof PotentialGroup) {
                AtomType[] types = ((PotentialGroup)link.potential).getAtomTypes(potential);
                if (types != null) {
                    return types;
                }
            }
        }
        return null;
    }
    
    /**
     * Returns an array containing the atom types for the molecules
     * corresponding to the given array of species.
     */
    protected AtomType[] moleculeTypes(Species[] species) {
        AtomType[] types = new AtomType[species.length];
        for(int i=0; i<species.length; i++) {
            types[i] = species[i].getFactory().getType();
        }
        return types;
    }
    
    private Species[] moleculeSpecies(AtomType[] types) {
        Species[] species = new Species[types.length];
        for (int i=0; i<types.length; i++) {
            species[i] = types[i].getSpecies();
        }
        return species;
    }
    
    protected void addPotential(Potential potential, AtomsetIteratorPDT iterator, AtomType[] types) {
        //the order of the given potential should be consistent with the order of the iterator
        if(potential.nBody() != iterator.nBody()) {
            throw new RuntimeException("Error: adding to PotentialGroup a potential and iterator that are incompatible");
        }
        //Set up to evaluate zero-body potentials last, since they may need other potentials
        //to be configured for calculation first
        if(potential instanceof Potential0) {//put zero-body potential at end of list
            if(last == null) {
                last = new PotentialLinker(potential, iterator, types, null);
                first = last;
            } else {
                last.next = new PotentialLinker(potential, iterator, types, null);
                last = last.next;
            }
        } else {//put other potentials at beginning of list
            first = new PotentialLinker(potential, iterator, types, first);
            if(last == null) last = first;
        }
        if (potential instanceof PotentialGroup) {
            ((PotentialGroup)potential).setPotentialMaster(this);
        }
    }

    /**
     * Removes given potential from the group.  No error is generated if
     * potential is not in group.
     */
    public synchronized void removePotential(Potential potential) {
        PotentialLinker previous = null;
        
        for(PotentialLinker link=first; link!=null; link=link.next) {
            if(link.potential == potential) {
                //found it
                if(previous == null) first = link.next;  //it's the first one
                else previous.next = link.next;          //it's not the first one
                //removing last; this works also if last was also first (then removing only, and set last to null)
                if(link == last) last = previous;
                return;
            }
            else if (link.potential instanceof PotentialGroup && 
                     ((PotentialGroup)link.potential).removePotential(potential)) {
                return;
            }
            previous = link;
        }
    }

 
    /**
     * @return Returns enabled flag.
     */
    public boolean isEnabled() {
        return enabled;
    }
    /**
     * Permits enabling/disabling of all potentials.  Default is enabled (true).
     * @param enabled flags if potentials are enabled.
     */
    public void setEnabled(boolean enabled) {
        this.enabled = enabled;
    }

    /**
     * Indicates that the specified potential should not contribute to potential
     * calculations. If potential is not in this group, no action is taken.
     */
    public void setEnabled(Potential potential, boolean enabled) {
        for(PotentialLinker link=first; link!=null; link=link.next) {
            if(link.potential == potential) {
                link.enabled = enabled;
                return;
            }
        }
    }
    
    /**
     * Returns true if the potential is in this group and has not been disabled
     * via a previous call to setEnabled; returns false otherwise.
     */
    public boolean isEnabled(Potential potential) {
        for(PotentialLinker link=first; link!=null; link=link.next) {
            if(link.potential == potential) {
                return link.enabled;
            }
        }
        return false;
    }
        
    /**
     * @return Returns the space.
     */
    public Space getSpace() {
        return space;
    }
    
    /**
     * Returns an array containing all molecular Potentials.
     */
    public Potential[] getPotentials() {
        int nPotentials=0;
        for(PotentialLinker link=first; link!=null; link=link.next) {
            nPotentials++;
        }
        Potential[] potentials = new Potential[nPotentials];
        int i=0;
        for(PotentialLinker link=first; link!=null; link=link.next) {
            potentials[i++] = link.potential;
        }
        return potentials;
    }
    
	protected PotentialMasterLrc lrcMaster;
	protected Phase mostRecentPhase = null;
	protected IteratorFactory iteratorFactory;

    protected PotentialLinker first, last;
    protected boolean enabled = true;
    protected final Space space;
    
    private static class AtomIterator0 extends AtomsetIteratorSinglet implements AtomsetIteratorPDT {
        AtomIterator0() {
            super(AtomSet.NULL);
        }
        public void setPhase(Phase phase) {}
        public void setTarget(Atom target) {}
        public void setDirection(IteratorDirective.Direction direction) {}
    }

    public static class PotentialLinker implements java.io.Serializable {
        public final Potential potential;
        public final AtomsetIteratorPDT iterator;
        public final AtomType[] types;
        public PotentialLinker next;
        public boolean enabled = true;
        //Constructors
        public PotentialLinker(Potential a, AtomsetIteratorPDT i, AtomType[] t, PotentialLinker l) {
            potential = a;
            iterator = i;
            next = l;
            if (t != null) {
                types = (AtomType[])t.clone();
            }
            else {
                types = null;
            }
        }
    }//end of PotentialLinker
}//end of PotentialMaster
    
