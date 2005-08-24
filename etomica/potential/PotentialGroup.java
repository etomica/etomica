package etomica.potential;

import java.util.Arrays;

import etomica.Phase;
import etomica.Space;
import etomica.Species;
import etomica.atom.AtomFilter;
import etomica.atom.AtomFilterTypeInstance;
import etomica.atom.AtomSet;
import etomica.atom.AtomType;
import etomica.atom.iterator.ApiBuilder;
import etomica.atom.iterator.AtomIteratorBasis;
import etomica.atom.iterator.AtomIteratorFiltered;
import etomica.atom.iterator.AtomsetIterator;
import etomica.atom.iterator.AtomsetIteratorBasisDependent;
import etomica.atom.iterator.AtomsetIteratorDirectable;
import etomica.atom.iterator.AtomsetIteratorSpeciesAgent;
import etomica.atom.iterator.IteratorDirective;
import etomica.atom.iterator.IteratorDirective.Direction;

/**
 * Collection of potentials that act between the atoms contained in
 * one or more groups of atoms.  This group iterates over all such atom-groups
 * assigned to it.  For each group it iterates over the potentials it contains,
 * instructing these sub-potentials to perform their calculations over the atoms
 * relevant to them in the groups.
 */
public class PotentialGroup extends Potential {
    
    /**
     * Makes a potential group defined on the position of nBody atom or atom groups.
     */
    public PotentialGroup(int nBody, Space space) {
        super(nBody, space);
    }

	/**
	 * Indicates if a given potential is a sub-potential of this group.
	 * @param potential the potential in question
	 * @return boolean true if potential has been added to this group
	 */
    public boolean contains(Potential potential) {
        for(PotentialLinker link=first; link!=null; link=link.next) {
            if(link.potential.equals(potential)) return true;
        }//end for
        return false;
    }
    
    /**
     * Adds the given potential and sets it up to apply to the atoms in the basis
     * having the given types.  Another addPotential method should be used if iteration
     * is not based on atom types. Length of types array must not exceed the order of
     * the potential (as given by the nBody method, and set in the constructor).
     * <br>
     * If length of types array is 1, this must be a one-body potential and iteration is
     * done over the atoms of the basis having the given type.  If types is of length 2
     * and this is a one-body potential, pairs are formed from atoms in the (single)
     * basis having the two given types (which might be the same); if this is a two-body
     * potential, pairs are formed from the first-type atoms taken from the first basis
     * atom, with the second-type atoms taken from the second basis.
     */
    public void addPotential(Potential potential, AtomType[] types, PotentialMaster potentialMaster) {
        if(this.nBody() > types.length) throw new IllegalArgumentException("Order of potential cannot exceed length of types array.");
        Arrays.sort(types);
        switch(types.length) {
            case 1:
                AtomFilter filter = new AtomFilterTypeInstance(types[0]);
                addPotential(potential, 
                        (AtomsetIteratorBasisDependent)AtomIteratorFiltered.makeIterator(new AtomIteratorBasis(),filter),types);
                break;
            case 2:
                if(this.nBody() == 1) {
                    addPotential(potential,
                            ApiBuilder.makeIntragroupTypeIterator(types),types);
                }
                else {//nBody == 2
                    addPotential(potential,
                            ApiBuilder.makeIntergroupTypeIterator(types),types);
                }
                break;
        }
        if(potential instanceof PotentialTruncated) {
            Potential0Lrc lrc = ((PotentialTruncated)potential).makeLrcPotential(types);
            if(lrc != null) {
                AtomsetIteratorSpeciesAgent iterator = new AtomsetIteratorSpeciesAgent(makeSpeciesArray(types));
                potentialMaster.lrcMaster().addPotential(lrc, iterator, null);
            }
        }
    }
 
    /**
     * Makes an array of species corresponding to the given array of types.
     */
    private Species[] makeSpeciesArray(AtomType[] types) {
        Species[] species = new Species[types.length];
        for(int i=0; i<types.length; i++) {
            species[i] = types[i].getSpecies();
        }
        return species;
    }
    
	/**
	 * Adds the given potential to this group, defining it to apply to the atoms
     * provided by the given basis-dependent iterator.  
	 */
	public synchronized void addPotential(Potential potential, AtomsetIteratorBasisDependent iterator) {
	    addPotential(potential,iterator,null);
    }
    
    private void addPotential(Potential potential, AtomsetIteratorBasisDependent iterator, AtomType[] types) {
        //the order of the given potential should be consistent with the order of the iterator
        if(potential.nBody() != iterator.nBody()) {
            throw new RuntimeException("Error: adding to PotentialGroup a potential and iterator that are incompatible");
        }
        //the given iterator should expect a basis of atoms equal in number to the order of this potential
        if(this.nBody() != iterator.basisSize()) {
            throw new RuntimeException("Error: adding an iterator that requires a basis size different from the nBody of the containing potential");
        }
        //put new potentials at beginning of list
        first = new PotentialLinker(potential, iterator, types, first);
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
	
	//TODO this needs some work
	public double energy(AtomSet basisAtoms) {
		if(basisAtoms.count() != this.nBody()) {
			throw new IllegalArgumentException("Error: number of atoms for energy calculation inconsistent with order of potential");
		}
		double sum = 0.0;
		for (PotentialLinker link=first; link!= null; link=link.next) {	
            if(!link.enabled) continue;
			//if(firstIterate) ((AtomsetIteratorBasisDependent)link.iterator).setDirective(id);
			link.iterator.setBasis(basisAtoms);
			link.iterator.reset();
			while(link.iterator.hasNext()) {
				sum += link.potential.energy(link.iterator.next());
			}
		}
		return sum;
	}
	
    /**
     * Returns the maximum of the range of all potentials in the group.
     */
	public double getRange() {
	    double range = 0;
        for(PotentialLinker link=first; link!=null; link=link.next) {
            if (link.potential.getRange() > range) {
                range = link.potential.getRange();
            }
        }
        return range;
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
				return;
			}//end if
			previous = link;
		}//end for
	}//end removePotential
    
    /**
     * Performs the specified calculation over the iterates given by the iterator,
     * using the directive to set up the iterators for the sub-potentials of this group.
     */
	//TODO consider what to do with sub-potentials after target atoms are reached
    public void calculate(AtomsetIterator iterator, IteratorDirective id, PotentialCalculation pc) {
	    	AtomSet targetAtoms = id.getTargetAtoms();
	    	IteratorDirective.Direction direction = id.direction();
			//loop over sub-potentials
	    	//TODO consider separate loops for targetable and directable
	 		for (PotentialLinker link=first; link!= null; link=link.next) {
				link.iterator.setTarget(targetAtoms);
	            if (link.iterator instanceof AtomsetIteratorDirectable) {
	                ((AtomsetIteratorDirectable)link.iterator).setDirection(direction);
	            }
			}
	    	iterator.reset();//loop over atom groups affected by this potential group
			while (iterator.hasNext()) {
	    		AtomSet basisAtoms = iterator.next();
	    		for (PotentialLinker link=first; link!= null; link=link.next) {
	                if(!link.enabled) continue;
	    			link.iterator.setBasis(basisAtoms);   			
	    			pc.doCalculation(link.iterator, id, link.potential);
	    		}
     	}
    }//end calculate
    
    public void setPhase(Phase phase) {
    	this.phase = phase;
  		for (PotentialLinker link=first; link!= null; link=link.next) {
			link.potential.setPhase(phase);
		}
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

	protected PotentialLinker first;
	protected Phase phase;

	protected static class PotentialLinker implements java.io.Serializable {
	    protected final Potential potential;
	    protected final AtomsetIteratorBasisDependent iterator;
        protected final AtomType[] types;
	    protected PotentialLinker next;
        protected boolean enabled = true;
	    //Constructors
	    public PotentialLinker(Potential a, AtomsetIteratorBasisDependent i, AtomType[] t, PotentialLinker l) {
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
	}

}//end PotentialGroup
    