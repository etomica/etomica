package etomica;


/**
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on Mar 3, 2005 by kofke
 */
public class AtomIndexManager {

    /**
     * 
     */
    public AtomIndexManager(int[] bitLength, int depth) {
        this.bitLength = (int[])bitLength.clone();
        cumulativeBitLength = calculateCumulativeBitLength(bitLength);
        bitShift = calculateBitShift(cumulativeBitLength);
        this.depth = depth;
        phaseMask = (power2(bitLength[1]) - 1) << (32-cumulativeBitLength[1]);
        speciesMask = (power2(bitLength[2]) - 1) << (32-cumulativeBitLength[2]);
        moleculeMask = (power2(bitLength[3]) - 1) << (32-cumulativeBitLength[3]);
        ordinalMask = (power2(bitLength[depth]) - 1) << (32-cumulativeBitLength[depth]);
    }
    
    private int power2(int n) {
        int power2 = 1;
        for(int i=0; i<n; i++) power2 *= 2;
        return power2;
    }
    // {speciesRoot, phases, species, molecules, groups, atoms}

    public AtomIndexManager makeChildManager() {
        return new AtomIndexManager(bitLength, depth+1);
    }
    
    public AtomIndexManager makeMoleculeIndexManager() {
        return new AtomIndexManager(bitLength, 3);
    }
    
    public int getDepth() {
        return depth;
    }
    
    public int shiftOrdinal(int ordinal) {
        if (((ordinal << bitShift[depth]) & ordinalMask) >>> bitShift[depth] != ordinal) {
            throw new RuntimeException(ordinal+" "+bitShift[depth]+" "+ordinalMask);
        }
        return ordinal << bitShift[depth];
    }
    
    public int getOrdinal(int atomIndex) {
        return (atomIndex & ordinalMask) >>> bitShift[depth];
    }
    
    public int getPhaseIndex(int unshiftedIndex) {
        return unshiftedIndex << bitShift[0];
    }
    
    public int getSpeciesIndex(int unshiftedIndex) {
        return unshiftedIndex << bitShift[1];
    }
    
    public int getMoleculeIndex(int unshiftedIndex) {
        return unshiftedIndex << bitShift[2];
    }
    
    public boolean samePhase(int index0, int index1) {
        return ((index0 ^ index1) & phaseMask) == 0;
    }

    public boolean sameSpecies(int index0, int index1) {
        return ((index0 ^ index1) & speciesMask) == 0;
    }

    public boolean sameMolecule(int index0, int index1) {
        return ((index0 ^ index1) & moleculeMask) == 0;
    }
    
    private static int[] calculateCumulativeBitLength(int[] array) {
        int[] newArray = (int[])array.clone();
        for(int i=1; i<newArray.length; i++) {
            newArray[i] += newArray[i-1];
        }
        if(newArray[newArray.length-1] != 31) {
            throw new RuntimeException("Improper setup of AtomTreeNode.MAX_ATOMS");
        }
        return newArray;
    }

    private static int[] calculateBitShift(int[] array) {
        int[] newArray = new int[array.length];
        for(int i=0; i<newArray.length; i++) {
            newArray[i] = 32 - array[i];
        }
        return newArray;
    }

    private final int[] bitLength;
    private final int[] cumulativeBitLength;
    private final int[] bitShift;
    private final int depth;
    private final int phaseMask;
    private final int speciesMask;
    private final int moleculeMask;
    private final int ordinalMask;
    
}
