/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.species;
import etomica.atom.IAtom;
import etomica.atom.IAtomType;
import etomica.api.IElement;
import etomica.api.IMolecule;
import etomica.atom.Atom;
import etomica.atom.AtomLeafDynamic;
import etomica.atom.AtomOriented;
import etomica.atom.AtomOrientedDynamic;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.IAtomTypeOriented;
import etomica.atom.Molecule;
import etomica.chem.elements.ElementSimple;
import etomica.config.ConformationLinear;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.util.Arrays;

/**
 * Species in which molecules are made of arbitrary number of atoms, with each
 * atom have a given atom type.  Atoms of a given type are grouped together.
 * The number of each type can be give explicitly or as a fraction of the
 * total.
 * 
 * @author David Kofke
 * @author Andrew Schultz
 */

public class SpeciesSpheresHetero extends Species {

    /**
     * Constructs instance with 0 components and total number of children 
     * equal to 1.  The actual atom types must be set before use.
     */
    public SpeciesSpheresHetero(Simulation sim, Space _space) {
        this(sim,_space, 0);
    }
    
    /**
     * Constructs instance with the given number of atom types.  Generic atom
     * types are created.
     */
    public SpeciesSpheresHetero(Simulation sim, Space _space, int nTypes) {
        this(_space, makeElements(sim,nTypes));
    }
    
    private static IElement[] makeElements(Simulation sim, int nTypes) {
        ElementSimple[] elements = new ElementSimple[nTypes];
        for (int i=0; i<elements.length; i++) {
            elements[i] = new ElementSimple(sim);
        }
        return elements;
    }
    
    /**
     * Constructs instance with the given elements.
     */
    public SpeciesSpheresHetero(Space _space, IElement[] leafElements) {
        this(_space, makeAtomTypes(leafElements));
    }
    
    protected static final AtomTypeLeaf[] makeAtomTypes(IElement[] leafElements) {
        AtomTypeLeaf[] types = new AtomTypeLeaf[leafElements.length];
        for (int i=0; i<types.length; i++) {
            types[i] = new AtomTypeLeaf(leafElements[i]);
        }
        return types;
    }
    
    /**
     * Constructs instance with the given atom types.
     */
    public SpeciesSpheresHetero(Space space, IAtomType[] atomTypes) {
        super();
        this.space = space;
        numberFraction = new double[atomTypes.length];
        childCount = new int[atomTypes.length];
        for (int i=0; i<atomTypes.length; i++) {
            addChildType(atomTypes[i]);
        }
        setConformation(new ConformationLinear(space));
    }
    
    public void setIsDynamic(boolean newIsDynamic) {
        isDynamic = newIsDynamic;
    }

    public boolean isDynamic() {
        return isDynamic;
    }

    /**
     * Constructs a new group containing a block of atoms for
     * each sub-type.
     */
    public IMolecule makeMolecule() {
        Molecule group = new Molecule(this, totalChildCount);
        // make block copolymers
        for (int i = 0; i < childTypes.length; i++) {
            for(int j = 0; j < childCount[i]; j++) {
                group.addChildAtom(makeLeafAtom(childTypes[i]));
            }
        }
        conformation.initializePositions(group.getChildList());
        return group;
    }

    protected IAtom makeLeafAtom(IAtomType leafType) {
        if (leafType instanceof IAtomTypeOriented) {
            return isDynamic ? new AtomOrientedDynamic(space, leafType)
                             : new AtomOriented(space, leafType);
        }
        return isDynamic ? new AtomLeafDynamic(space, leafType)
                         : new Atom(space, leafType);
    }

    public int getNumComponents() {
        return childTypes.length;
    }

    /**
     * Sets the number fraction of each child atom type made by
     * this factory.  The number of child atoms for each type
     * is calculated from the number fraction for that type and
     * the total number of child atoms.  The number fractions normalized and
     * the number of children is calculated in order of increasing number 
     * fraction.
     * 
     * @throws IllegalArgumentException
     *             if childFactory.length != childCount.length
     *             if any number fraction is <0 or >1
     *             if all number fractions=0
     */
    public void setNumberFraction(double[] newNumberFraction) {
        if (newNumberFraction.length != childTypes.length) {
            throw new IllegalArgumentException("the number of numberFracions must be equal to the number of child factories");
        }
        numberFraction = newNumberFraction.clone();
        childCount = new int[numberFraction.length];
        
        double totalFraction = 0;
        //uninitialize childCount, determine totalFraction so that 
        //numberFraction can be renormalized
        for (int i=0; i<childCount.length; i++) {
            childCount[i] = -1;
            if (numberFraction[i] < 0.0 || numberFraction[i] > 1.0) {
                throw new IllegalArgumentException("number fraction must be between 0 and 1 (inclusive)");
            }
            totalFraction += numberFraction[i];
        }
        if (totalFraction == 0) {
            throw new IllegalArgumentException("All number fractions cannot be 0!");
        }
        
        // determine the number of each child type, starting with the type with
        // the smallest number fraction.  Calculate number for each type based 
        // on the fraction remaining and the number of atoms remaining so that
        // the number of children is equal to totalChildCount
        int atomsRemaining = totalChildCount;
        for (int i=0; i<childTypes.length; i++) {
            int k = -1;
            double minFraction = Double.POSITIVE_INFINITY;
            for (int j=0; j<childTypes.length; j++) {
                if (childCount[j] != -1) {
                    // already did this type
                    continue;
                }
                if (numberFraction[j] < minFraction) {
                    k = j;
                    minFraction = numberFraction[j];
                }
            }
            childCount[k] = (int)Math.round(atomsRemaining*(numberFraction[k]/totalFraction));
            totalFraction -= numberFraction[k];
            atomsRemaining -= childCount[k];
        }
    }
    
    /**
     * Returns the number fraction of each child type created by this factory.
     */
    public double[] getNumberFraction() {
        return numberFraction.clone();
    }
    
    /**
     * Directly sets the number of child atoms of each type created by this 
     * factory. Invoking this method also sets numberFraction and 
     * totalChildCount based on the given parameters.
     * 
     * @throws IllegalArgumentException
     *             if childFactory.length != childCount.length
     *             if any childCount is < 0
     */
    public void setChildCount(int[] newChildCount) {
        if (childTypes.length != newChildCount.length) {
            throw new IllegalArgumentException("Number of child factories must equal length of newChildCount.  " +
                    "Call setChildFactory first");
        }
        childCount = newChildCount.clone();
        totalChildCount = 0;
        for (int i=0; i<childCount.length; i++) {
            if (childCount[i] < 0) {
                throw new IllegalArgumentException("Each child count must be positive");
            }
            totalChildCount += childCount[i];
        }
        
        //calculate number fraction from childCount
        for (int i=0; i<childCount.length; i++) {
            numberFraction[i] = (double)childCount[i]/totalChildCount;
        }
    }

    /**
     * Returns the number of child atoms of each type created by this factory.
     */
    public int[] getChildCount() {
        return childCount.clone();
    }
    
    /**
     * Sets the factories that make the child atoms of this factory's atom.
     * If the number of factories changes, the number fractions are set so 
     * that there is an equal amount of each child.  The caller is responsible
     * for ensuring that the AtomTypes for the child factories are children
     * of this AtomFactory's AtomType.
     *
     * @throws IllegalArgumentException
     *             if newChildFactory is an empty array
     */
    public void setChildTypes(IAtomType[] newchildTypes) {
        for (int i=0; i<childTypes.length; i++) {
            removeChildType(childTypes[i]);
        }
        for (int i=0; i<childTypes.length; i++) {
            addChildType(newchildTypes[i]);
        }
        if (numberFraction.length != childTypes.length) {
            double[] fraction = new double[childTypes.length];
            java.util.Arrays.fill(fraction, 1.0/childTypes.length);
            setNumberFraction(fraction);
        }
    }

    /**
     * Adds the given factory as a child of this factory.  The caller is 
     * responsible for ensuring that the AtomType for the child factory is a
     * child of this AtomFactory's AtomType.
     */
    public void addChildType(IAtomType newLeafType) {
        super.addChildType(newLeafType);
        if (childTypes.length > 1) {
            // assume fraction = 0 for new childFactory
            numberFraction = Arrays.resizeArray(numberFraction,numberFraction.length+1);
            childCount = Arrays.resizeArray(childCount,childCount.length+1);
        }
        else {
            setNumberFraction(new double[]{1.0});
        }
    }
    
    /**
     * Removes the given factory as a child of this factory.  The caller is 
     * responsible for removing the AtomType for the child factory from its
     * parent AtomType.
     */
    public void removeChildType(IAtomType oldLeafType) {
        super.removeChildType(oldLeafType);
        if (childTypes.length > 0) {
            double[] newNumberFraction = new double[numberFraction.length-1];
            System.arraycopy(numberFraction,0,newNumberFraction,0,index);
            System.arraycopy(numberFraction,index+1,newNumberFraction,index,numberFraction.length-index-1);
            // setNumberFraction will renormalize the fractions
            setNumberFraction(newNumberFraction);
        }
        else {
            numberFraction = new double[0];
            childCount = new int[0];
        }
    }
    
    /**
     * Sets the total number of child atoms created by this factory.  The 
     * number of each atom type is recalculated based on the previous
     * number fractions.
     *
     * @throws IllegalArgumentException
     *             if numChildren < 0
     */
    public void setTotalChildren(int numChildren) {
        if (numChildren < 0) {
            throw new IllegalArgumentException("Number of children must not be negative");
        }
        totalChildCount = numChildren;
        if (numberFraction.length > 0) {
            setNumberFraction(numberFraction);
        }
    }

    public int getNumLeafAtoms() {
        int total = 0;
        for (int i=0; i<childCount.length; i++) {
            total += childCount[i];
        }
        return total;
    }

    private static final long serialVersionUID = 1L;
    protected Space space;
    protected boolean isDynamic;
    protected double[] numberFraction;
    protected int[] childCount;
    protected int totalChildCount;
}
