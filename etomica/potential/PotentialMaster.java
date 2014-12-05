/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.api.IPotentialAtomic;
import etomica.api.IPotentialMaster;
import etomica.api.IPotentialMolecular;
import etomica.api.ISpecies;
import etomica.atom.AtomArrayList;
import etomica.atom.MoleculeArrayList;
import etomica.atom.iterator.AtomsetIteratorBasisDependent;
import etomica.atom.iterator.IteratorDirective;
import etomica.atom.iterator.IteratorFactory;
import etomica.atom.iterator.MoleculeIteratorAll;
import etomica.atom.iterator.MoleculesetIteratorPDT;
import etomica.chem.models.Model;
import etomica.chem.models.Model.PotentialAndIterator;


/**
 * Manager of all potentials in simulation.
 * Most calls to compute the energy or other potential calculations begin
 * with the calculate method of this class.  It then passes the calculation 
 * on to the contained potentials.
 *
 * @author David Kofke
 */
public class PotentialMaster implements IPotentialMaster {
    
    public PotentialMaster() {
        this(IteratorFactory.INSTANCE);
    } 
    
    public PotentialMaster(IteratorFactory iteratorFactory) {
        this.iteratorFactory = iteratorFactory;
        potentialList = new ArrayList<PotentialLinker>();
    }
    
    /* (non-Javadoc)
	 * @see etomica.potential.IPotentialMaster#lrcMaster()
	 */
	 public PotentialMasterLrc lrcMaster() {
		if(lrcMaster == null) lrcMaster = new PotentialMasterLrc();
		return lrcMaster;
	 }

     /* (non-Javadoc)
	 * @see etomica.potential.IPotentialMaster#makePotentialGroup(int)
	 */
     public PotentialGroup makePotentialGroup(int nBody) {
         return new PotentialGroup(nBody);
     }

     /* (non-Javadoc)
	 * @see etomica.potential.IPotentialMaster#calculate(etomica.box.Box, etomica.atom.iterator.IteratorDirective, etomica.potential.PotentialCalculation)
	 */
    public void calculate(IBox box, IteratorDirective id, PotentialCalculation pc) {
        if(!enabled) return;
    	IMolecule targetMolecule = id.getTargetMolecule();
    	IAtom targetAtomLeaf = id.getTargetAtom();
    	if (targetAtomLeaf != null) {
    	    targetMolecule = targetAtomLeaf.getParentGroup();
    	}

        for(PotentialLinker link : potentialList) {
    	    if(!link.enabled) continue;
    	    final MoleculesetIteratorPDT atomIterator = link.iterator;
    	    final IPotentialMolecular potential = link.potential;
	        atomIterator.setBox(box);
	        potential.setBox(box);
    	    atomIterator.setTarget(targetMolecule);
    	    atomIterator.setDirection(id.direction());
    	    if (pc instanceof PotentialCalculationMolecular) {
    	        atomIterator.reset();
                for (IMoleculeList atoms = atomIterator.next(); atoms != null;
                     atoms = atomIterator.next()) {
                    ((PotentialCalculationMolecular)pc).doCalculation(atoms, potential);
                }
    	    }
    	    else if (potential instanceof PotentialGroup) {
    	        ((PotentialGroup)potential).calculate(atomIterator, id.direction(), targetAtomLeaf, pc);
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
        if (getPotential(new ISpecies[]{newModel.getSpecies()}) != null) {
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
    public void addPotential(IPotentialMolecular potential, ISpecies[] species) {
    	if (potential.nBody() == 0) {
    		addPotential(potential, new MoleculeIterator0(),null);
    	}
        else if (potential.nBody() == Integer.MAX_VALUE) {
            addPotential(potential, new MoleculeIteratorAll(species), null);
        }
    	else if (potential.nBody() != species.length) {
    		throw new IllegalArgumentException("Illegal species length");
    	}
        else {
            MoleculesetIteratorPDT iterator = iteratorFactory.makeMoleculeIterator(species);
            addPotential(potential, iterator, species);
        }
    }
    
    /* (non-Javadoc)
	 * @see etomica.potential.IPotentialMaster#addPotential(etomica.api.IPotential, etomica.api.IAtomType[])
	 */
    public void addPotential(IPotentialAtomic potential, IAtomType[] atomTypes) {
        if (potential.nBody() != Integer.MAX_VALUE && potential.nBody() != atomTypes.length) {
            throw new IllegalArgumentException("nBody of potential must match number of atom types");
        }
        
    	Arrays.sort(atomTypes);
        // depth of molecules
        ISpecies[] parentAtomTypes = new ISpecies[atomTypes.length];
        for (int i=0; i<atomTypes.length; i++) {
            parentAtomTypes[i] = atomTypes[i].getSpecies();
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
    public void potentialAddedNotify(IPotentialAtomic subPotential, PotentialGroup pGroup) {
        // do nothing.  this is here for subclasses to override
    }
    
    /* (non-Javadoc)
	 * @see etomica.potential.IPotentialMaster#getPotential(etomica.api.IAtomType[])
	 */
    public PotentialGroup getPotential(ISpecies[] types) {
        for(PotentialLinker link : potentialList) {
            if (link.potential instanceof PotentialGroup) {
                if(Arrays.equals(types,link.types)) {
                    return (PotentialGroup)link.potential;
                }
            }
        }
        return null;
    }
    
    public ISpecies[] getSpecies(IPotentialMolecular potential) {
        for(PotentialLinker link : potentialList) {
            if (link.potential == potential) {
                return link.types;
            }
        }
        return null;
    }

    public void addPotential(IPotentialMolecular potential, MoleculesetIteratorPDT iterator, ISpecies[] types) {
        //the order of the given potential should be consistent with the order of the iterator
        if(potential.nBody() != iterator.nBody()) {
            throw new RuntimeException("Error: adding to PotentialGroup a potential and iterator that are incompatible");
        }
        //Set up to evaluate zero-body potentials last, since they may need other potentials
        //to be configured for calculation first
        potentialList.add(new PotentialLinker(potential, iterator, types));
        if (potential instanceof PotentialGroup) {
            ((PotentialGroup)potential).setPotentialMaster(this);
        }
    }

    /* (non-Javadoc)
	 * @see etomica.potential.IPotentialMaster#removePotential(etomica.api.IPotential)
	 */
    public synchronized void removePotential(IPotentialMolecular potential) {
        
        for(PotentialLinker link : potentialList) {
            if(link.potential == potential) {
                potentialList.remove(link);
                return;
            }
        }
        throw new RuntimeException("potential not found");
    }
    
    public synchronized void removePotential(IPotentialAtomic potential) {
        
        for(PotentialLinker link : potentialList) {
            if (link.potential instanceof PotentialGroup && ((PotentialGroup)link.potential).contains(potential)) {
                ((PotentialGroup)link.potential).removePotential(potential);
                return;
            }
        }
        throw new RuntimeException("potential not found");
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
    public void setEnabled(IPotentialMolecular potential, boolean enabled) {
        for(PotentialLinker link : potentialList) {
            if(link.potential == potential) {
                link.enabled = enabled;
                return;
            }
        }
    }
    
    public void setEnabled(IPotentialAtomic potential, boolean enabled) {
        for(PotentialLinker link : potentialList) {
            if (link.potential instanceof PotentialGroup && ((PotentialGroup)link.potential).contains(potential)) {
                ((PotentialGroup)link.potential).setEnabled(potential, enabled);
                return;
            }
        }
    }
    
    public boolean isEnabled(IPotentialAtomic potential) {
        for(PotentialLinker link : potentialList) {
            if (link.potential instanceof PotentialGroup && ((PotentialGroup)link.potential).contains(potential)) {
                return ((PotentialGroup)link.potential).isEnabled(potential);
            }
        }
        throw new RuntimeException("I don't know about "+potential);
    }
    
    /* (non-Javadoc)
	 * @see etomica.potential.IPotentialMaster#isEnabled(etomica.potential.Potential)
	 */
    public boolean isEnabled(IPotentialMolecular potential) {
        for(PotentialLinker link : potentialList) {
            if(link.potential == potential) {
                return link.enabled;
            }
        }
        throw new RuntimeException("I don't know about "+potential);
    }
    
    /* (non-Javadoc)
	 * @see etomica.potential.IPotentialMaster#getPotentials()
	 */
    public IPotentialMolecular[] getPotentials() {
        IPotentialMolecular[] potentials = new IPotentialMolecular[potentialList.size()];
        int i=0;
        for(PotentialLinker link : potentialList) {
            potentials[i++] = link.potential;
        }
        return potentials;
    }
    
	protected PotentialMasterLrc lrcMaster;
	protected IteratorFactory iteratorFactory;

    protected List<PotentialLinker> potentialList;
    protected boolean enabled = true;

    public static class MoleculeIterator0 implements MoleculesetIteratorPDT {
        public MoleculeIterator0() {}
        public void setBox(IBox box) {}
        public void setTarget(IMolecule target) {}
        public void setDirection(IteratorDirective.Direction direction) {}
        public int nBody() {return 0;}
        public IMoleculeList next() {
            if (finished) {
                return null;
            }
            finished = true;
            return list;
        }
        public void reset() {finished = false;}
        public int size() {return 1;}
        public void unset() {finished = true;}
        protected final MoleculeArrayList list = new MoleculeArrayList(0);
        protected boolean finished;
    }
    //Added
    public static class AtomIterator0 implements AtomsetIteratorBasisDependent {
        public AtomIterator0() {}
        public void setTarget(IAtom target) {}
        public int nBody() {return 0;}
        public IAtomList next() {
            if (finished) {
                return null;
            }
            finished = true;
            return list;
        }
        public void reset() {finished = false;}
        public int size() {return 1;}
        public void unset() {finished = true;}
        protected final AtomArrayList list = new AtomArrayList(0);
        protected boolean finished;
		public void setBasis(IMoleculeList atoms) {}
		public int basisSize() {return 0;}
		public boolean haveTarget(IAtom target) {return false;}
    }


    public static class PotentialLinker {
        public final IPotentialMolecular potential;
        public final MoleculesetIteratorPDT iterator;
        public final ISpecies[] types;
        public boolean enabled = true;
        //Constructors
        public PotentialLinker(IPotentialMolecular a, MoleculesetIteratorPDT i, ISpecies[] t) {
            potential = a;
            iterator = i;
            if (t != null) {
                types = t.clone();
            }
            else {
                types = null;
            }
        }
    }//end of PotentialLinker
}//end of PotentialMaster
    
