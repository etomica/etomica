package etomica.chem.models;

import etomica.api.ISimulation;
import etomica.api.ISpecies;
import etomica.atom.iterator.ApiBuilder;
import etomica.chem.elements.Element;
import etomica.chem.elements.ElementSimple;
import etomica.config.Conformation;
import etomica.potential.Potential2;
import etomica.space.Space;
import etomica.species.SpeciesSpheres;

/**
 * Model class for a simple chain molecule (homopolymer)
 * @author Andrew Schultz
 */
public class ModelChain extends Model {

    /**
     * Constructs a default ModelChain instance.  The number of atoms per chain
     * must be set before making the Species.  The conformation and element can
     * also be set if desired.
     */
    public ModelChain(Space _space) {
        speciesMade = false;
        space = _space;
    }

    /**
     * Set the bonding potential for adjacent atoms in the chain molecule.
     * The given potential must not be range-dependent.  This method must be
     * invoked before making the Species and must not be invoked after creating
     * the Species.
     */
    public void setBondingPotential(Potential2 newBondingPotential) {
        if (speciesMade) {
            throw new RuntimeException("Species already created");
        }
        if (newBondingPotential.getRange() < Double.POSITIVE_INFINITY) {
            throw new IllegalArgumentException("bonding potential must not be range-dependent");
        }
        bondingPotential = newBondingPotential;
    }
    
    /**
     * Returns the bonding potential used for adjacent atoms in the chain
     * molecule.
     */
    public Potential2 getBondingPotential() {
        return bondingPotential;
    }
    
    /**
     * Sets the number of atoms per chain molecule.  The number of atoms must
     * be positive.  This method may not be invoked after creating the Species.
     */
    public void setNumAtoms(int newNumAtoms) {
        if (speciesMade) {
            throw new RuntimeException("Species already created");
        }
        if (newNumAtoms < 1) {
            throw new IllegalArgumentException("Number of atoms must be positive");
        }
        numAtoms = newNumAtoms;
    }
    
    /**
     * Returns the number of atoms per chain molecule.
     */
    public int getNumAtoms() {
        return numAtoms;
    }

    /**
     * Sets the Conformation for the molecule.  If no Conformation
     * is set, a ConformationLinear will be used.
     */
    public void setConformation(Conformation newConformation) {
        if (speciesMade) {
            throw new RuntimeException("Species already created");
        }
        conformation = newConformation;
    }
    
    /**
     * Returns the Conformation for the molecule, or null if no Conformation
     * has been set (a ConformationLinear will be used).
     */
    public Conformation getConformation() {
        return conformation;
    }

    /**
     * Sets the Element used for the leaf atoms in the chain.
     */
    public void setElement(Element newElement) {
        if (speciesMade) {
            throw new RuntimeException("Species already created");
        }
        element = newElement;
    }
    
    /**
     * Returns the Element used for the leaf atoms in the chain, or null if no
     * Element has been set (an ElementSimple will be used).
     */
    public Element getElement() {
        return element;
    }
    
    protected void initPotentials(ISimulation sim) {
        // we already have our bonding potential, so do nothing
    }

    public PotentialAndIterator[] getPotentials() {
        return new PotentialAndIterator[]{new PotentialAndIterator(
                bondingPotential,ApiBuilder.makeAdjacentPairIterator())};
    }

    protected ISpecies makeSpeciesInternal(ISimulation sim) {
        if (bondingPotential == null) {
            throw new RuntimeException("Please set the bonding potential before" +
                    " creating the Species");
        }
        
        if (numAtoms == 0) {
            throw new RuntimeException("You must first set the number of leaf atoms");
        }
        
        if (element == null) {
            setElement(new ElementSimple(sim));
        }
        
        if (conformation == null) {
            ISpecies species = new SpeciesSpheres(sim, space, numAtoms, element);
            setConformation(species.getMoleculeType().getConformation());
        }
        
        speciesMade = true;
        return new SpeciesSpheres(sim, numAtoms, element, conformation, space);
    }

    private static final long serialVersionUID = 1L;
    protected Potential2 bondingPotential;
    protected int numAtoms;
    protected Conformation conformation;
    protected Element element;
    protected boolean speciesMade;
    private final Space space;
}
