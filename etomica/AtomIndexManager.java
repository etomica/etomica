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
        int rootMask = 1 << 31;
        phaseIndexMask = ((power2(bitLength[1]) - 1) << bitShift[1]);
        speciesIndexMask = ((power2(bitLength[2]) - 1) << bitShift[2]);
        moleculeIndexMask = ((power2(bitLength[3]) - 1) << bitShift[3]); 
        ordinalMask = (power2(bitLength[depth]) - 1) << bitShift[depth];
        samePhaseMask = phaseIndexMask | rootMask;
        sameSpeciesMask = speciesIndexMask | rootMask;
        sameMoleculeMask = moleculeIndexMask | rootMask | samePhaseMask | sameSpeciesMask;
//        System.out.println("depth, bitLength: "+depth+Arrays.toString(bitLength));
//        System.out.println("samePhaseMask: "+Integer.toBinaryString(samePhaseMask));
//        System.out.println("sameSpeciesMask: "+Integer.toBinaryString(sameSpeciesMask));
//        System.out.println("sameMoleculeMask: "+Integer.toBinaryString(sameMoleculeMask));
//        System.out.println("ordinalMask: "+Integer.toBinaryString(ordinalMask|rootMask));
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
        return (unshiftedIndex & phaseIndexMask) >>> bitShift[1];
    }
    
    public int getSpeciesIndex(int unshiftedIndex) {
        return (unshiftedIndex & speciesIndexMask) >>> bitShift[2];
    }
    
    public int getMoleculeIndex(int unshiftedIndex) {
        return (unshiftedIndex & moleculeIndexMask) >>> bitShift[3];
    }
    
    public boolean samePhase(int index0, int index1) {
        return ((index0 ^ index1) & samePhaseMask) == 0;
    }

    public boolean sameSpecies(int index0, int index1) {
        return ((index0 ^ index1) & sameSpeciesMask) == 0;
    }

    public boolean sameMolecule(int index0, int index1) {
        return ((index0 ^ index1) & sameMoleculeMask) == 0;
    }
    
    private static int[] calculateCumulativeBitLength(int[] array) {
        int[] newArray = (int[])array.clone();
        for(int i=1; i<newArray.length; i++) {
            newArray[i] += newArray[i-1];
        }
        if(newArray[newArray.length-1] != 32) {
            throw new RuntimeException("Improper setup of Bitlength");
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
    private final int samePhaseMask;
    private final int sameSpeciesMask;
    private final int sameMoleculeMask;
    private final int phaseIndexMask;
    private final int speciesIndexMask;
    private final int moleculeIndexMask;
    private final int ordinalMask;
    
}
