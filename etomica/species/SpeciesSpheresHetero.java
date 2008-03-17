package etomica.species;
import java.lang.reflect.Constructor;

import etomica.api.IMolecule;
import etomica.api.ISimulation;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomLeafDynamic;
import etomica.atom.AtomPositionGeometricCenter;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.AtomTypeMolecule;
import etomica.atom.AtomTypeSphere;
import etomica.atom.IAtomLeaf;
import etomica.atom.Molecule;
import etomica.chem.elements.Element;
import etomica.chem.elements.ElementSimple;
import etomica.space.Space;
import etomica.util.Arrays;

/**
 * Species in which molecules are made of arbitrary number of spheres,
 * with each sphere having the same mass and size (same type).
 * 
 * @author David Kofke
 */

public class SpeciesSpheresHetero extends Species {

    /**
     * Constructs instance with 0 components and total number of children 
     * equal to 1.  The actual number of components must be set in the factory
     * (AtomFactoryHetero) before use.  The actual number of desired children 
     * can also be set in the factory.
     */
    public SpeciesSpheresHetero(ISimulation sim, Space _space) {
        this(sim,_space, 0);
    }
    
    /**
     * Constructs instance with the given number of components and 
     * total number of children equal to 1.  The actual number of desired 
     * desired children can be set in the factory (AtomFactoryHetero) after
     * construction.
     */
    public SpeciesSpheresHetero(ISimulation sim, Space _space, int nComponents) {
        this(sim, _space, makeElements(sim,nComponents));
    }
    
    private static Element[] makeElements(ISimulation sim, int nComponents) {
        ElementSimple[] elements = new ElementSimple[nComponents];
        for (int i=0; i<elements.length; i++) {
            elements[i] = new ElementSimple(sim);
        }
        return elements;
    }
    
    /**
     * Constructs instance with the given leaf elements and 
     * total number of children equal to 1.  The actual number of desired 
     * desired children can be set in the factory (AtomFactoryHetero) after
     * construction.
     */
    public SpeciesSpheresHetero(ISimulation sim, Space _space, Element[] leafElements) {
        this(_space, sim.isDynamic(), makeAtomTypeSpheres(leafElements));
    }
    
    protected static final AtomTypeSphere[] makeAtomTypeSpheres(Element[] leafElements) {
        AtomTypeSphere[] types = new AtomTypeSphere[leafElements.length];
        for (int i=0; i<types.length; i++) {
            types[i] = new AtomTypeSphere(leafElements[i]);
        }
        return types;
    }
    
    public SpeciesSpheresHetero(Space space, boolean isDynamic, AtomTypeLeaf[] atomTypes) {
        super(new AtomTypeMolecule(new AtomPositionGeometricCenter(space)));
        this.space = space;
        this.isDynamic = isDynamic;
        this.leafTypes = (AtomTypeLeaf[])atomTypes.clone();
        numberFraction = new double[atomTypes.length];
        childCount = new int[atomTypes.length];
        for (int i=0; i<atomTypes.length; i++) {
            atomType.addChildType(atomTypes[i]);
        }
    }
    
    /**
     * Constructs a new group containing a block of atoms for
     * each sub-type.
     */
    public IMolecule makeMolecule() {
        isMutable = false;
        Molecule group = new Molecule(atomType);
        // make block copolymers
        for (int i = 0; i < leafTypes.length; i++) {
            for(int j = 0; j < childCount[i]; j++) {
                group.addChildAtom(makeLeafAtom(leafTypes[i]));
            }
        }
        return group;
    }

    protected IAtomLeaf makeLeafAtom(AtomTypeLeaf leafType) {
        return isDynamic ? new AtomLeafDynamic(space, leafType)
                         : new AtomLeaf(space, leafType);
    }

    public int getNumComponents() {
        return leafTypes.length;
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
        if (!isMutable) {
            throw new IllegalStateException("Factory is not mutable");
        }
        if (newNumberFraction.length != leafTypes.length) {
            throw new IllegalArgumentException("the number of numberFracions must be equal to the number of child factories");
        }
        numberFraction = (double[])newNumberFraction.clone();
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
        for (int i=0; i<leafTypes.length; i++) {
            int k = -1;
            double minFraction = Double.POSITIVE_INFINITY;
            for (int j=0; j<leafTypes.length; j++) {
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
        if (!isMutable) {
            throw new IllegalStateException("Factory is not mutable");
        }
        if (leafTypes.length != newChildCount.length) {
            throw new IllegalArgumentException("Number of child factories must equal length of newChildCount.  " +
                    "Call setChildFactory first");
        }
        childCount = (int[])newChildCount.clone();
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
    public void setLeafTypes(AtomTypeLeaf[] newLeafTypes) {
        if (!isMutable) {
            throw new IllegalStateException("Factory is not mutable");
        }
        for (int i=0; i<leafTypes.length; i++) {
            atomType.removeChildType(leafTypes[i]);
        }
        leafTypes = (AtomTypeLeaf[])newLeafTypes.clone();
        for (int i=0; i<leafTypes.length; i++) {
            atomType.addChildType(newLeafTypes[i]);
        }
        if (numberFraction.length != leafTypes.length) {
            double[] fraction = new double[leafTypes.length];
            java.util.Arrays.fill(fraction, 1.0/leafTypes.length);
            setNumberFraction(fraction);
        }
    }

    /**
     * Adds the given factory as a child of this factory.  The caller is 
     * responsible for ensuring that the AtomType for the child factory is a
     * child of this AtomFactory's AtomType.
     */
    public void addLeafType(AtomTypeLeaf newLeafType) {
        if (!isMutable) {
            throw new IllegalStateException("Factory is not mutable");
        }
        atomType.addChildType(newLeafType);
        leafTypes = (AtomTypeLeaf[])Arrays.addObject(leafTypes,newLeafType);
        if (leafTypes.length > 1) {
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
    public boolean removeLeafType(AtomTypeLeaf oldLeafType) {
        if (!isMutable) {
            throw new IllegalStateException("Factory is not mutable");
        }
        int typeIndex = -1;
        for (int i=0; i<leafTypes.length; i++) {
            if (leafTypes[i] == oldLeafType) {
                typeIndex = i;
                break;
            }
        }
        if (typeIndex == -1) {
            return false;
        }
        leafTypes = (AtomTypeLeaf[])Arrays.removeObject(leafTypes, oldLeafType);
        if (leafTypes.length > 0) {
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
        return true;
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
        if (!isMutable) {
            throw new IllegalStateException("Factory is not mutable");
        }
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

    
    public SpeciesSignature getSpeciesSignature() {
        Constructor constructor = null;
        try {
            constructor = this.getClass().getConstructor(new Class[]{ISimulation.class,Element[].class});
        }
        catch(NoSuchMethodException e) {
            System.err.println("you have no constructor.  be afraid");
        }
        Element[] leafElements = new Element[leafTypes.length];
        for (int i=0; i<leafTypes.length; i++) {
            leafElements[i] = leafTypes[i].getElement();
        }
        //XXX broken
        return new SpeciesSignature(constructor,new Object[]{leafElements});
    }
    
    private static final long serialVersionUID = 1L;
    protected Space space;
    protected boolean isDynamic;
    protected double[] numberFraction;
    protected int[] childCount;
    protected int totalChildCount;
    protected AtomTypeLeaf[] leafTypes;
}
