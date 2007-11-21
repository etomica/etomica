package etomica.atom;

import etomica.config.Conformation;
import etomica.config.ConformationLinear;
import etomica.simulation.ISimulation;
import etomica.space.Space;
import etomica.util.Arrays;

/**
 * Builds an atom group that comprises a set of differently-formed atoms or
 * atomgroups. Each child atom is constructed by a different atom factory. 
 * Position definition is the center-of-mass.
 * 
 * @author David Kofke
 */

public class AtomFactoryHetero extends AtomFactory {

    public AtomFactoryHetero(ISimulation sim) {
        this(sim.getSpace(), new ConformationLinear(sim));
    }

    public AtomFactoryHetero(Space space, Conformation conformation) {
        super(new AtomTypeGroup(new AtomPositionCOM(space)));
        ((AtomTypeGroup)atomType).setConformation(conformation);
        childFactory = new AtomFactory[0];
        numberFraction = new double[0];
        childCount = new int[0];
    }

    /**
     * Constructs a new group containing a block of atoms for
     * each sub-type.
     */
    public IAtom makeAtom() {
        isMutable = false;
        AtomGroup group = new AtomGroup(atomType);
        // make block copolymers
        for (int i = 0; i < childFactory.length; i++) {
            for(int j = 0; j < childCount[i]; j++) {
                group.addChildAtom(childFactory[i].makeAtom());
            }
        }
        return group;
    }

    public int getNumComponents() {
        return childFactory.length;
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
        if (newNumberFraction.length != childFactory.length) {
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
        for (int i=0; i<childFactory.length; i++) {
            int k = -1;
            double minFraction = Double.POSITIVE_INFINITY;
            for (int j=0; j<childFactory.length; j++) {
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
        return (double[])numberFraction.clone();
        // don't use numberFraction since it won't necessarily be equal to 
        // the actual fraction generated for each type
//        double[] fraction = new double[childCount.length];
//        for (int i=0; i<fraction.length; i++) {
//            fraction[i] = (double)childCount[i]/totalChildCount;
//        }
//        return fraction;
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
        if (childFactory.length != newChildCount.length) {
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
        return (int[])childCount.clone();
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
    public void setChildFactory(AtomFactory[] newChildFactory) {
        if (!isMutable) {
            throw new IllegalStateException("Factory is not mutable");
        }
        childFactory = (AtomFactory[])newChildFactory.clone();
        if (numberFraction.length != childFactory.length) {
            double[] fraction = new double[childFactory.length];
            java.util.Arrays.fill(fraction, 1.0/childFactory.length);
            setNumberFraction(fraction);
        }
    }

    /**
     * Adds the given factory as a child of this factory.  The caller is 
     * responsible for ensuring that the AtomType for the child factory is a
     * child of this AtomFactory's AtomType.
     */
    public void addChildFactory(AtomFactory newChildFactory) {
        if (!isMutable) {
            throw new IllegalStateException("Factory is not mutable");
        }
        childFactory = (AtomFactory[])Arrays.addObject(childFactory,newChildFactory);
        if (childFactory.length > 1) {
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
    public boolean removeChildFactory(AtomFactory newChildFactory) {
        if (!isMutable) {
            throw new IllegalStateException("Factory is not mutable");
        }
        int index = -1;
        for (int i=0; i<childFactory.length; i++) {
            if (childFactory[i] == newChildFactory) {
                index = i;
                break;
            }
        }
        if (index == -1) {
            return false;
        }
        childFactory = (AtomFactory[])Arrays.removeObject(childFactory,newChildFactory);
        if (childFactory.length > 0) {
            double[] newNumberFraction = new double[numberFraction.length-1];
            System.arraycopy(numberFraction,0,newNumberFraction,0,index);
            System.arraycopy(numberFraction,index+1,newNumberFraction,index,numberFraction.length-index-1);
            setNumberFraction(newNumberFraction);
        }
        else {
            numberFraction = new double[0];
            childCount = new int[0];
        }
        return true;
    }
    
    /**
     * Returns the array of subfactories that produces each of the 
     * atoms in the group made by this factory.
     */
    public AtomFactory[] getChildFactory() {
        return (AtomFactory[])childFactory.clone();
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

    public int getNumChildAtoms() {
        return totalChildCount;
    }
    
    public int getNumTreeAtoms() {
        int total = 1;
        for (int i=0; i<childCount.length; i++) {
            total += childCount[i] * childFactory[i].getNumTreeAtoms();
        }
        return total;
    }

    public int getNumLeafAtoms() {
        int total = 0;
        for (int i=0; i<childCount.length; i++) {
            total += childCount[i] * childFactory[i].getNumLeafAtoms();
        }
        return total;
    }

    private static final long serialVersionUID = 1L;
    protected AtomFactory[] childFactory;
    protected int[] childCount;
    protected int totalChildCount;
    protected double[] numberFraction;
}
