package etomica.potential;

import java.util.Arrays;

import etomica.api.IAtom;
import etomica.api.IAtomSet;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.api.IPotential;
import etomica.api.IPotentialMaster;
import etomica.api.ISpecies;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.AtomsetArray;
import etomica.atom.iterator.AtomIteratorAll;
import etomica.atom.iterator.AtomsetIteratorPDT;
import etomica.atom.iterator.AtomsetIteratorSinglet;
import etomica.atom.iterator.IteratorDirective;
import etomica.atom.iterator.IteratorFactory;
import etomica.chem.models.Model;
import etomica.chem.models.Model.PotentialAndIterator;
import etomica.space.Space;


/**
 * Manager of all potentials in simulation.
 * Most calls to compute the energy or other potential calculations begin
 * with the calculate method of this class.  It then passes the calculation 
 * on to the contained potentials.
 *
 * @author David Kofke
 */
public class PotentialMaster implements java.io.Serializable, IPotentialMaster {
    
    public PotentialMaster(Space space) {
        this(space,IteratorFactory.INSTANCE);
    } 
    
    public PotentialMaster(Space space, IteratorFactory iteratorFactory) {
        this.space = space;
        this.iteratorFactory = iteratorFactory;
    }
    
    /* (non-Javadoc)
	 * @see etomica.potential.IPotentialMaster#lrcMaster()
	 */
	 public PotentialMasterLrc lrcMaster() {
		if(lrcMaster == null) lrcMaster = new PotentialMasterLrc(space);
		return lrcMaster;
	 }

     /* (non-Javadoc)
	 * @see etomica.potential.IPotentialMaster#makePotentialGroup(int)
	 */
     public PotentialGroup makePotentialGroup(int nBody) {
         return new PotentialGroup(nBody,space);
     }

     /* (non-Javadoc)
	 * @see etomica.potential.IPotentialMaster#calculate(etomica.box.Box, etomica.atom.iterator.IteratorDirective, etomica.potential.PotentialCalculation)
	 */
    public void calculate(IBox box, IteratorDirective id, PotentialCalculation pc) {
        if(!enabled) return;
    	IAtom targetAtom = id.getTargetAtom();
    	mostRecentBox = box;

        for(PotentialLinker link=first; link!=null; link=link.next) {
    	    if(!link.enabled) continue;
    	    final AtomsetIteratorPDT atomIterator = link.iterator;
    	    final IPotential potential = link.potential;
	        atomIterator.setBox(box);
	        potential.setBox(box);
    	    atomIterator.setTarget(targetAtom);
    	    atomIterator.setDirection(id.direction());
    	    if (potential instanceof PotentialGroup) {
    	        ((PotentialGroup)potential).calculate(atomIterator, id, pc);
    	    }
    	    else {
    	        atomIterator.reset();
                for (IAtomSet atoms = atomIterator.next(); atoms != null;
                     atoms = atomIterator.next()) {
                    pc.doCalculation(atoms, potential);
                }
    	    }
        }
        
        if(lrcMaster != null) {
            lrcMaster.calculate(box, id, pc);
        }
    }
    
    /* (non-Javadoc)
	 * @see etomica.potential.IPotentialMaster#addModel(etomica.chem.models.Model)
	 */
    public void addModel(Model newModel) {
        if (getPotential(new IAtomType[]{newModel.getSpecies().getMoleculeType()}) != null) {
            throw new IllegalArgumentException(newModel+" has already been added");
        }
        PotentialAndIterator[] potentialsAndIterators = newModel.getPotentials();
        PotentialGroup pGroup = makePotentialGroup(1);
        for (int i=0; i<potentialsAndIterators.length; i++) {
            pGroup.addPotential(potentialsAndIterators[i].getPotential(),
                    potentialsAndIterators[i].getIterator());
        }
        addPotential(pGroup, new ISpecies[]{newModel.getSpecies()});
    }
    
    /* (non-Javadoc)
	 * @see etomica.potential.IPotentialMaster#addPotential(etomica.api.IPotential, etomica.species.ISpecies[])
	 */
    public void addPotential(IPotential potential, ISpecies[] species) {
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
    
    /* (non-Javadoc)
	 * @see etomica.potential.IPotentialMaster#addPotential(etomica.api.IPotential, etomica.api.IAtomType[])
	 */
    public void addPotential(IPotential potential, IAtomType[] atomTypes) {
    	if (potential.nBody() == 0) {
    		addPotential(potential, new AtomIterator0(),null);
    		return;
    	}    	
        else if (potential.nBody() != Integer.MAX_VALUE && potential.nBody() != atomTypes.length) {
            throw new IllegalArgumentException("nBody of potential must match number of atom types");
        }
        
    	Arrays.sort(atomTypes);
        // depth of molecules
        boolean haveLeafTypes = false;
        for (int i=0; i<atomTypes.length; i++) {
            if (atomTypes[i] instanceof AtomTypeLeaf) {
                haveLeafTypes = true;
            }
        }
        if (!haveLeafTypes) {
            addPotential(potential,moleculeSpecies(atomTypes));
            return;
        }
        IAtomType[] parentAtomTypes = new IAtomType[atomTypes.length];
        for (int i=0; i<atomTypes.length; i++) {
            if (atomTypes[i] instanceof AtomTypeLeaf) {
                parentAtomTypes[i] = ((AtomTypeLeaf)atomTypes[i]).getParentType();
            }
            else {
                parentAtomTypes[i] = atomTypes[i];
            }
        }
        // look for a PotentialGroup that applies to parentAtomTypes
        PotentialGroup pGroup = getPotential(parentAtomTypes);
        if (pGroup == null) { // didn't find an appropriate potentialgroup
            pGroup = makePotentialGroup(potential.nBody());
            addPotential(pGroup,parentAtomTypes);
        }
        pGroup.addPotential(potential,atomTypes);
    }
    
    /* (non-Javadoc)
	 * @see etomica.potential.IPotentialMaster#potentialAddedNotify(etomica.api.IPotential, etomica.potential.PotentialGroup)
	 */
    public void potentialAddedNotify(IPotential subPotential, PotentialGroup pGroup) {
        // do nothing.  this is here for subclasses to override
    }
    
    /* (non-Javadoc)
	 * @see etomica.potential.IPotentialMaster#getPotential(etomica.api.IAtomType[])
	 */
    public PotentialGroup getPotential(IAtomType[] types) {
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
    
    /* (non-Javadoc)
	 * @see etomica.potential.IPotentialMaster#getAtomTypes(etomica.api.IPotential)
	 */
    public IAtomType[] getAtomTypes(IPotential potential) {
        for(PotentialLinker link=first; link!=null; link=link.next) {
            if (link.potential == potential) {
                return link.types;
            }
            if (link.potential instanceof PotentialGroup) {
                IAtomType[] types = ((PotentialGroup)link.potential).getAtomTypes(potential);
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
    protected IAtomType[] moleculeTypes(ISpecies[] species) {
        IAtomType[] types = new IAtomType[species.length];
        for(int i=0; i<species.length; i++) {
            types[i] = species[i].getMoleculeType();
        }
        return types;
    }
    
    private ISpecies[] moleculeSpecies(IAtomType[] types) {
        ISpecies[] species = new ISpecies[types.length];
        for (int i=0; i<types.length; i++) {
            species[i] = types[i].getSpecies();
        }
        return species;
    }
    
    public void addPotential(IPotential potential, AtomsetIteratorPDT iterator, IAtomType[] types) {
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

    /* (non-Javadoc)
	 * @see etomica.potential.IPotentialMaster#removePotential(etomica.api.IPotential)
	 */
    public synchronized void removePotential(IPotential potential) {
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

 
    /* (non-Javadoc)
	 * @see etomica.potential.IPotentialMaster#isEnabled()
	 */
    public boolean isEnabled() {
        return enabled;
    }
    /* (non-Javadoc)
	 * @see etomica.potential.IPotentialMaster#setEnabled(boolean)
	 */
    public void setEnabled(boolean enabled) {
        this.enabled = enabled;
    }

    /* (non-Javadoc)
	 * @see etomica.potential.IPotentialMaster#setEnabled(etomica.potential.Potential, boolean)
	 */
    public void setEnabled(IPotential potential, boolean enabled) {
        for(PotentialLinker link=first; link!=null; link=link.next) {
            if(link.potential == potential) {
                link.enabled = enabled;
                return;
            }
        }
    }
    
    /* (non-Javadoc)
	 * @see etomica.potential.IPotentialMaster#isEnabled(etomica.potential.Potential)
	 */
    public boolean isEnabled(IPotential potential) {
        for(PotentialLinker link=first; link!=null; link=link.next) {
            if(link.potential == potential) {
                return link.enabled;
            }
        }
        return false;
    }
        
    /* (non-Javadoc)
	 * @see etomica.potential.IPotentialMaster#getSpace()
	 */
    public Space getSpace() {
        return space;
    }
    
    /* (non-Javadoc)
	 * @see etomica.potential.IPotentialMaster#getPotentials()
	 */
    public IPotential[] getPotentials() {
        int nPotentials=0;
        for(PotentialLinker link=first; link!=null; link=link.next) {
            nPotentials++;
        }
        IPotential[] potentials = new Potential[nPotentials];
        int i=0;
        for(PotentialLinker link=first; link!=null; link=link.next) {
            potentials[i++] = link.potential;
        }
        return potentials;
    }
    
    private static final long serialVersionUID = 1L;
	protected PotentialMasterLrc lrcMaster;
	protected IBox mostRecentBox = null;
	protected IteratorFactory iteratorFactory;

    protected PotentialLinker first, last;
    protected boolean enabled = true;
    protected final Space space;

    public static class AtomIterator0 extends AtomsetIteratorSinglet implements AtomsetIteratorPDT {
        private static final long serialVersionUID = 1L;
        public AtomIterator0() {
            super(new AtomsetArray(0));
        }
        public void setBox(IBox box) {}
        public void setTarget(IAtom target) {}
        public void setDirection(IteratorDirective.Direction direction) {}
    }

    public static class PotentialLinker implements java.io.Serializable {
        private static final long serialVersionUID = 1L;
        public final IPotential potential;
        public final AtomsetIteratorPDT iterator;
        public final IAtomType[] types;
        public PotentialLinker next;
        public boolean enabled = true;
        //Constructors
        public PotentialLinker(IPotential a, AtomsetIteratorPDT i, IAtomType[] t, PotentialLinker l) {
            potential = a;
            iterator = i;
            next = l;
            if (t != null) {
                types = (IAtomType[])t.clone();
            }
            else {
                types = null;
            }
        }
    }//end of PotentialLinker
}//end of PotentialMaster
    
