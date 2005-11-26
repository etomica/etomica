package etomica.atom;

import java.util.Arrays;

import etomica.config.Conformation;
import etomica.config.ConformationLinear;
import etomica.data.DataSourceCOM;
import etomica.simulation.Simulation;
import etomica.space.CoordinateFactoryNull;
import etomica.space.Space;
import etomica.species.Species;

/**
 * Builds an atom group that comprises a set of differently-formed atoms or
 * atomgroups. Each child atom is constructed by a different atom factory. 
 * Position definition is the center-of-mass.
 * 
 * @author David Kofke
 */

public class AtomFactoryHetero extends AtomFactory {

    public AtomFactoryHetero(Simulation sim,
            AtomSequencerFactory sequencerFactory, AtomTypeGroup parentType) {
        this(sim.space, sequencerFactory, parentType, new ConformationLinear(sim));
    }

    public AtomFactoryHetero(Space space,
            AtomSequencerFactory sequencerFactory, AtomTypeGroup parentType,
            Conformation config) {
        this(space, sequencerFactory, parentType, AtomTreeNodeGroup.FACTORY,
                config);
    }

    public AtomFactoryHetero(Space space,
            AtomSequencerFactory sequencerFactory, AtomTypeGroup parentType,
            AtomTreeNodeFactory nodeFactory, Conformation config) {
        super(new CoordinateFactoryNull(), new AtomTypeGroup(parentType, new DataSourceCOM(space)),
                sequencerFactory, nodeFactory);
        conformation = config;
        childFactory = new AtomFactory[0];
        numberFraction = new double[0];
        childCount = new int[0];
    }

    public void setSpecies(Species species) {
        atomType.setSpecies(species);
        for (int i = 0; i < childFactory.length; i++) {
            childFactory[i].setSpecies(species);
        }
    }

    /**
     * Constructs a new group containing a block of atoms for
     * each sub-type.
     */
    public Atom makeAtom() {
        Atom group = newParentAtom();
        AtomTreeNodeGroup node = (AtomTreeNodeGroup) group.node;
        // make block copolymers
        for (int i = 0; i < childFactory.length; i++) {
            for(int j = 0; j < childCount[i]; j++) {
                Atom childAtom = childFactory[i].makeAtom();
                childAtom.node.setParent(node);
            }
        }
        return group;
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
            double minFraction = 2.0;
            for (int j=0; j<childFactory.length; j++) {
                if (childCount[i] != -1) {
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
        // don't use numberFraction since it won't necessarily be equal to 
        // the actual fraction generated for each type
        double[] fraction = new double[childCount.length];
        for (int i=0; i<fraction.length; i++) {
            fraction[i] = (double)childCount[i]/totalChildCount;
        }
        return fraction;
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
     * that there is an equal amount of each child.
     *
     * @throws IllegalArgumentException
     *             if newChildFactory is an empty array
     */
    public void setChildFactory(AtomFactory[] newChildFactory) {
        if (childFactory.length != 0) {
            throw new IllegalStateException("You can set the child factory only once!");
        }
        childFactory = (AtomFactory[])newChildFactory.clone();
        if (numberFraction.length != childFactory.length) {
            double[] fraction = new double[childFactory.length];
            Arrays.fill(fraction, 1.0/childFactory.length);
            setNumberFraction(fraction);
        }
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
        if (numChildren < 0) {
            throw new IllegalArgumentException("Number of children must not be negative");
        }
        totalChildCount = numChildren;
        if (numberFraction.length > 0) {
            setNumberFraction(numberFraction);
        }
    }

    /**
     * Returns the total number of child atoms created by this factory.
     */
    public int getTotalChildren() {
        return totalChildCount;
    }

    private AtomFactory[] childFactory;
    private int[] childCount;
    private int totalChildCount;
    private double[] numberFraction;
}
