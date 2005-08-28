package etomica.atom;


/**
 * Performs manipulations and interrogations related to the index assigned to
 * each atom. The atom index is held in the node field of the atom. It is an
 * integer treated as a set of bits, with different bit segments coding
 * information about where the atom is located in the atom tree. The methods of
 * this class provide a convenient means for obtaining information from this
 * code. Instances of this class are held in the AtomType instance referenced by
 * the type field of each atom.
 * <p> 
 * In addition, the AtomIndexManager holds an index that establishes a parallel
 * hierarchy to the atom hierarchy.  The index is coded in the same way as that
 * for the atoms (although the ordinal values will be smaller, as there are
 * many atom instances with the same AtomIndexManager).  This index is accessed
 * via the getTypeIndex method.
 * <p>
 * Instances of this class are held as a field of AtomType.
 * 
 * @see AtomTreeNode
 */

/*
 * History Created on Mar 3, 2005 by kofke
 */
public class AtomIndexManager implements java.io.Serializable {

    /**
     * @param bitLength
     *            Indicates how the atom index is coded, by specifying the
     *            number of bits used to locate an atom at each depth of the
     *            tree. This array is specified to the constructor of the
     *            governing Simulation instance, and it is propagated down to
     *            the index manager via this constructor. The bitLength array is
     *            not changed after the Simulation is constructed. The first
     *            element of the array should be 1, and the value of this bit
     *            indicates whether the atom is in the tree or not (0 means it
     *            is not, 1 means it is, so if the index is read as an integer
     *            it will be negative if the atom is in the tree hierarchy);
     *            subsequent bits indicate the species master (phase), species
     *            agent, molecule, then groups/atoms in the molecule as
     *            appropriate to the system being simulated.
     * @param depth
     *            specifies the depth of this atom in the atom tree hierarchy.
     *            The depth and the bitLength array together specify the set of
     *            bits in the atom index that code for the ordinal index of the
     *            atom. The ordinal index is a simple integer that is assigned
     *            sequentially (1, 2, 3, etc) to each atom as it is added to its
     *            parent atom. Each atom has an ordinal index, and this with the
     *            ordinal indexes of all parents up the atom tree, form the
     *            atom's index that is managed by this class.
     *  
     */
    private AtomIndexManager(int[] bitLength, int depth, int parentIndex, int ordinal) {
        this.bitLength = (int[]) bitLength.clone();
        cumulativeBitLength = calculateCumulativeBitLength(bitLength);
        bitShift = calculateBitShift(cumulativeBitLength);
        this.depth = depth;
        int rootMask = 1 << 31;
        phaseIndexMask = ((power2(bitLength[1]) - 1) << bitShift[1]);
        speciesIndexMask = ((power2(bitLength[2]) - 1) << bitShift[2]);
        moleculeIndexMask = ((power2(bitLength[3]) - 1) << bitShift[3]);
        ordinalMask = (power2(bitLength[depth]) - 1) << bitShift[depth];
        indexMask = (power2(cumulativeBitLength[depth])-1) << bitShift[depth];//indexMask has 1's in for all bits significant to the index
        samePhaseMask = phaseIndexMask | rootMask;
        sameSpeciesMask = speciesIndexMask | rootMask;
        sameMoleculeMask = moleculeIndexMask | rootMask | samePhaseMask
                | sameSpeciesMask;
        typeIndex = parentIndex + shiftOrdinal(ordinal);
//        System.out.println(Integer.toBinaryString(typeIndex));
//        System.out.println("depth, bitLength,cumulativeBitLength,bitShift: "+depth+Arrays.toString(bitLength)+Arrays.toString(cumulativeBitLength)+Arrays.toString(bitShift));
//        System.out.println("samePhaseMask: "+Integer.toBinaryString(samePhaseMask));
//        System.out.println("sameSpeciesMask: "+Integer.toBinaryString(sameSpeciesMask));
//        System.out.println("sameMoleculeMask: "+Integer.toBinaryString(sameMoleculeMask));
//        System.out.println("ordinalMask: "+Integer.toBinaryString(ordinalMask|rootMask));
//        System.out.println("  indexMask: "+Integer.toBinaryString(indexMask));
    }

    //convenience method; return 2^n
    private int power2(int n) {
        return 1 << n;
    }

    // {speciesRoot, phases, species, molecules, groups, atoms}
    /*
     * SpeciesRoot has unique AtomType
     * All SpeciesMasters in a simulation share the same AtomType (root.childType)
     * Each Species has its own indexManager, which it gives to all its SpeciesAgents
     * All molecules of a particular species have the same indexManager
     */

    /**
     * Returns an AtomIndexManager instance that would be used by a child of
     * this manager's atom. Simply constructs and returns an index manager with
     * the same bitLength array as this, with its depth field incremented by 1,
     * and with an ordinal index that is increment with each call to this method.
     * Method exists to permit construction of appropriate index manger having
     * same bitLength array as this, without exposing the bitLength array.
     */
    AtomIndexManager makeChildManager() {
        return new AtomIndexManager(bitLength, depth + 1, typeIndex, ++childCount);
    }

    /**
     * Constructs an index manager for the SpeciesRoot AtomType.  Called in
     * SpeciesRoot constructor.
     */
    static AtomIndexManager makeRootIndexManager(int[] bitLength) {
        return new AtomIndexManager(bitLength, 0, 0, 1);
    }
    
    /**
     * Special-use method to make index manager. Constructs an index manager 
     * at the molecule depth, with typeIndex = 0, and ordinal of 1.  Called by
     * AtomType constructor if parent type is null.
     */
    public static AtomIndexManager makeSimpleIndexManager(int[] bitLength) {
        return new AtomIndexManager(bitLength, 3, 0, 1);
    }

    /**
     * @return the value of the depth field of this manager.
     */
    public int getDepth() {
        return depth;
    }

    /**
     * Bit-shifts the given integer by the amount needed to put it in the bit
     * location appropriate to the depth of the atom. The atom's index would
     * then be given by adding the returned value to the index of the atom's
     * parent. The ordinal can be recovered from the index via the getOrdinal
     * method.
     */
    public int shiftOrdinal(int ordinal) {
        if (((ordinal << bitShift[depth]) & ordinalMask) >>> bitShift[depth] != ordinal) {
            throw new RuntimeException(ordinal + " " + bitShift[depth] + " "
                    + ordinalMask);
        }
        return ordinal << bitShift[depth];
    }

    /**
     * Returns the atom's ordinal index by decoding its atomIndex. Reverses the
     * action of shiftOrdinal.
     */
    public int getOrdinal(int atomIndex) {
        return (atomIndex & ordinalMask) >>> bitShift[depth];
    }

    /**
     * Decodes an atom's index to determine the index of the phase it is in.
     */
    public int getPhaseIndex(int unshiftedIndex) {
        return (unshiftedIndex & phaseIndexMask) >>> bitShift[1];
    }

    /**
     * Decodes an atom's index to determine the index of the species it is part
     * of.
     */
    public int getSpeciesIndex(int unshiftedIndex) {
        return (unshiftedIndex & speciesIndexMask) >>> bitShift[2];
    }

    /**
     * Decodes an atom's index to determine the ordinal index of the molecule it
     * is part of.
     */
    public int getMoleculeIndex(int unshiftedIndex) {
        return (unshiftedIndex & moleculeIndexMask) >>> bitShift[3];
    }

    /**
     * Returns true if the given indices correspond to atoms that are in the 
     * same phase.
     * @param index0 index of an atom, as given by the index() method of its node.
     * @param index1 index of another atom, as given by the index() method of its node.
     * @return true if atoms are in the same phase (or are the same atom).
     */
    public boolean samePhase(int index0, int index1) {
        return ((index0 ^ index1) & samePhaseMask) == 0;
    }

    /**
     * Returns true if the given indices correspond to atoms that are part of the 
     * same species.
     * @param index0 index of an atom, as given by the index() method of its node.
     * @param index1 index of another atom, as given by the index() method of its node.
     * @return true if atoms are in the same species (or are the same atom).
     */
    public boolean sameSpecies(int index0, int index1) {
        return ((index0 ^ index1) & sameSpeciesMask) == 0;
    }

    /**
     * Returns true if the given indices correspond to atoms that are part of the 
     * same molecule.
     * @param index0 index of an atom, as given by the index() method of its node.
     * @param index1 index of another atom, as given by the index() method of its node.
     * @return true if atoms are in the same molecule (or are the same atom).
     */
    public boolean sameMolecule(int index0, int index1) {
        return ((index0 ^ index1) & sameMoleculeMask) == 0;
    }
    
    /**
     * Returns true if the indices are for atoms that have the
     * same ancestry to the depth of this manager.  If true, and one of the 
     * indices is for an atom at this depth, then the other index is for
     * an atom descended from it. 
     * @param index0
     * @param index1
     * @return
     */
    public boolean sameAncestry(int index0, int index1) {
        return ((index0 ^ index1) & indexMask) == 0;
    }
    
    /**
     * Returns true if an atom of this manager's type is descended from
     * an atom (any atom) having the type of given manager.
     */
    public boolean isDescendedFrom(AtomIndexManager anotherManager) {
        return anotherManager.sameAncestry(typeIndex, anotherManager.typeIndex);
    }

    /**
     * Returns the index associated with this manager, and which is used to
     * locate this manager's type relative to the other atom types.  This index
     * is set at construction.
     */
    public int getTypeIndex() {
        return typeIndex;
    }

    //convenience method used by constructor
    private static int[] calculateCumulativeBitLength(int[] array) {
        int[] newArray = (int[]) array.clone();
        for (int i = 1; i < newArray.length; i++) {
            newArray[i] += newArray[i - 1];
        }
        if (newArray[newArray.length - 1] != 32) {
            throw new RuntimeException("Improper setup of Bitlength");
        }
        return newArray;
    }

    //convenience method used by constructor
    private static int[] calculateBitShift(int[] array) {
        int[] newArray = new int[array.length];
        for (int i = 0; i < newArray.length; i++) {
            newArray[i] = 32 - array[i];
        }
        return newArray;
    }
    
    private final int[] bitLength;
    private final int[] cumulativeBitLength;
    private final int[] bitShift;
    private final int depth;
    private final int samePhaseMask;
    private final int sameSpeciesMask;
    private final int sameMoleculeMask;
    private final int phaseIndexMask;
    private final int speciesIndexMask;
    private final int moleculeIndexMask;
    private final int ordinalMask;
    private final int indexMask;

    private final int typeIndex;
    private int childCount = 0;
}
