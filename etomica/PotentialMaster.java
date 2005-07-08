package etomica;

import java.util.Arrays;

import etomica.atom.AtomSequencerFactory;
import etomica.atom.iterator.AtomsetIteratorMolecule;
import etomica.atom.iterator.AtomsetIteratorSinglet;
import etomica.potential.Potential0;
import etomica.potential.Potential0Lrc;
import etomica.potential.PotentialCalculation;
import etomica.potential.PotentialTruncated;


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
        this(space,IteratorFactorySimple.INSTANCE);
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

     /**
      * Performs the given PotentialCalculation on the atoms of the given Phase.
      * Sets the phase for all molecule iterators and potentials, sets target
      * and direction for iterators as specified by given IteratorDirective,
      * and applies doCalculation of given PotentialCalculation with the iterators
      * and potentials.
      */
    public void calculate(Phase phase, IteratorDirective id, PotentialCalculation pc) {
    	if(!enabled) return;
    	AtomSet targetAtoms = id.getTargetAtoms();
    	boolean phaseChanged = (phase != mostRecentPhase);
    	mostRecentPhase = phase;
    	for(PotentialLinker link=first; link!=null; link=link.next) {
            if(!link.enabled) continue;
			if(phaseChanged) {
				link.iterator.setPhase(phase);
				link.potential.setPhase(phase);
			}
			link.iterator.setTarget(targetAtoms);
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
    public void setSpecies(Potential potential, Species[] species) {
    	if (potential.nBody() == 0) {
    		addPotential(potential, new AtomIterator0(),null);
    	}
    	else if (potential.nBody() != species.length) {
    		throw new IllegalArgumentException("Illegal species length");
    	}
        else {
            AtomsetIteratorMolecule iterator = iteratorFactory.makeMoleculeIterator(species);
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
     * another.  Potential group heirarchy will be constructed as needed above
     * the level of the given atom types.
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
            depth[i] = atomTypes[i].getIndexManager().getDepth();
            if (depth[i] > maxDepth) maxDepth = depth[i];
        }
        if (maxDepth == 3) {
            setSpecies(potential,moleculeSpecies(atomTypes));
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
            pGroup = new PotentialGroup(atomTypes.length,space);
            addPotential(pGroup,parentAtomTypes);
        }
        pGroup.addPotential(potential,atomTypes,this);
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
    
    protected void addPotential(Potential potential, AtomsetIteratorMolecule iterator, AtomType[] types) {
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
    }

    /**
     * Removes given potential from the group.  No error is generated if
     * potential is not in group.
     */
    public synchronized void removePotential(Potential potential) {
        PotentialLinker previous = null;
        for(PotentialLinker link=first; link!=null; link=link.next) {
            if(link.potential == potential) {//found it
                if(previous == null) first = link.next;  //it's the first one
                else previous.next = link.next;          //it's not the first one
                if(link == last) last = previous; //removing last; this works also if last was also first (then removing only, and set last to null)
                return;
            }//end if
            previous = link;
        }//end for
    }//end removePotential

 
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
        
    public AtomSequencerFactory sequencerFactory() {
        return iteratorFactory.moleculeSequencerFactory();
    }

    
    /**
     * @return Returns the space.
     */
    public Space getSpace() {
        return space;
    }
    
	protected PotentialMasterLrc lrcMaster;
	protected Phase mostRecentPhase = null;
	protected IteratorFactory iteratorFactory;

    protected PotentialLinker first, last;
    protected boolean enabled = true;
    protected final Space space;
    
    private static class AtomIterator0 extends AtomsetIteratorSinglet implements AtomsetIteratorMolecule {
        AtomIterator0() {
            super(AtomSet.NULL);
        }
        public void setPhase(Phase phase) {}
        public void setTarget(AtomSet target) {}
        public void setDirection(IteratorDirective.Direction direction) {}
    }

    public static class PotentialLinker implements java.io.Serializable {
        public final Potential potential;
        public final AtomsetIteratorMolecule iterator;
        public final AtomType[] types;
        public PotentialLinker next;
        public boolean enabled = true;
        //Constructors
        public PotentialLinker(Potential a, AtomsetIteratorMolecule i, AtomType[] t, PotentialLinker l) {
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
    
