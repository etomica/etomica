/*
 * History
 * Created on Dec 3, 2004 by kofke
 */
package etomica;

/**
 * Provides static methods for construction of some useful pair iterators.
 */
public final class ApiBuilder {

    /**
     * Private constructor to prevent instantiation.
     */
    private ApiBuilder() {}

    /**
     * Pair iterator that constructs pairs from the childlist atoms of a
     * basis, each with the atom adjacent to it in the list. If no target is
     * set, returns each such pair once.  If a target is set, returns the target
     * atom with the atom uplist from it (if direction is not set, or is set to UP) or down
     * from it (if direction is set to DOWN), or both (if direction is set to null).
     */
    public static ApiIntragroup makeAdjacentPairIterator() {
        AtomIteratorSequencerBonded aiInner = new AtomIteratorSequencerBonded();
        aiInner.setNumToSkip(1);
        return new ApiIntragroup(new ApiInnerVariable(
                                        new AtomIteratorBasis(),aiInner));
    }

    /**
     * Pair iterator that constructs pairs from the childlist atoms of a
     * basis, each with all atoms that are not adjacent to it in the list. If no target is
     * set, returns each such pair once.  If a target is set, returns the target
     * atom with the non-adjacent atoms uplist from it (if direction is not set, 
     * or is set to UP) or down from it (if direction is set to DOWN), or both 
     * (if direction is set to null).
     */
    public static ApiIntragroup makeNonAdjacentPairIterator() {
        AtomIteratorSequencerList aiInner = new AtomIteratorSequencerList();
        aiInner.setNumToSkip(2);
        return new ApiIntragroup(new ApiInnerVariable(
                                        new AtomIteratorBasis(),aiInner));
    }
    
    /**
     * Makes a pair iterator that forms pairs from the atoms of two different lists.
     * To set the list, access the inner and outer iterators, thus:<br>
     * <code>
     * ApiInnerFixed pairIterator = makeInterlistIterator();
     * ((AtomIteratorListSimple)pairIterator.getInnerIterator).setList(innerList);
     * ((AtomIteratorListSimple)pairIterator.getOuterIterator).setList(outerList);
     * </code>
     */
    public static ApiInnerFixed makeInterlistIterator() {
        return new ApiInnerFixed(new AtomIteratorListSimple(), new AtomIteratorListSimple());
    }
}
