/*
 * History
 * Created on Dec 3, 2004 by kofke
 */
package etomica.atom.iterator;

import etomica.atom.AtomFilter;
import etomica.atom.AtomFilterTypeInstance;
import etomica.atom.AtomType;

/**
 * Provides static methods for construction of some useful pair iterators.
 */
public final class ApiBuilder implements java.io.Serializable {

    /**
     * Private constructor to prevent instantiation.
     */
    private ApiBuilder() {
    }

    /**
     * Pair iterator that constructs pairs from the childlist atoms of a basis,
     * each with the atom adjacent to it in the list. If no target is set,
     * returns each such pair once. If a target is set, returns the target atom
     * with the atom uplist from it (if direction is not set, or is set to UP)
     * or down from it (if direction is set to DOWN), or both (if direction is
     * set to null).
     */
    public static ApiIntragroup makeAdjacentPairIterator() {
        AtomIteratorSequenceAdjacent aiInner = new AtomIteratorSequenceAdjacent();
        return new ApiIntragroup(new ApiInnerVariable(new AtomIteratorBasis(),
                aiInner));
    }

    /**
     * Pair iterator that constructs pairs from the childlist atoms of a basis,
     * each with all atoms that are not adjacent to it in the list. If no target
     * is set, returns each such pair once. If a target is set, returns the
     * target atom with the non-adjacent atoms uplist from it (if direction is
     * not set, or is set to UP) or down from it (if direction is set to DOWN),
     * or both (if direction is set to null).
     */
    public static ApiIntragroup makeNonAdjacentPairIterator() {
        AtomIteratorSequenceDirectable aiInner = new AtomIteratorSequenceDirectable();
        aiInner.setNumToSkip(2);
        return new ApiIntragroup(new ApiInnerVariable(new AtomIteratorBasis(),
                aiInner));
    }

    /**
     * Returns an intergroup iterator that filters the iterates so that only
     * those having the given type instances are returned. The given types are
     * applied corresponding to the pair of basis atoms identifed when the
     * iterator's setBasis method is invoked. Child atoms of the first basis
     * atom only having type = types[0] are given, in pairs with child atoms of
     * the second basis atom only having type = types[1].
     * 
     * @throws IllegalArgumentException if the types array is not of length 2
     */
    public static ApiIntergroup makeIntergroupTypeIterator(AtomType[] types) {
        if (types.length != 2)
            throw new IllegalArgumentException(
                    "Incorrect number of types; must be 2");
        AtomFilter typeFilter0 = new AtomFilterTypeInstance(types[0]);
        AtomFilter typeFilter1 = new AtomFilterTypeInstance(types[1]);
        AtomIterator outer = AtomIteratorFiltered.makeIterator(
                new AtomIteratorBasis(), typeFilter0);
        AtomIterator inner = AtomIteratorFiltered.makeIterator(
                new AtomIteratorBasis(), typeFilter1);
        return new ApiIntergroup(new ApiInnerFixed(outer, inner));
    }

    /**
     * Returns an intragroup iterator that filters the iterates so that only
     * those having the given type instances are returned. The given types are
     * applied correspond to the pair of basis atoms identifed when the
     * iterator's setBasis method is invoked. Child atoms of the first basis
     * atom only having type = types[0] are given, in pairs with child atoms of
     * the second basis atom only having type = types[1]. If the two given types
     * are the different instances, an intergroup iterator is returned. An
     * exception is thrown if the types array is not of length 2.
     */
    public static AtomsetIteratorBasisDependent makeIntragroupTypeIterator(
            AtomType[] types) {
        if (types.length != 2)
            throw new IllegalArgumentException(
                    "Incorrect number of types; must be 2");
        if (types[0] == types[1]) {
            AtomFilter typeFilter = new AtomFilterTypeInstance(types[0]);
            AtomIterator outer = AtomIteratorFiltered.makeIterator(
                    new AtomIteratorBasis(), typeFilter);
            AtomIteratorFiltered inner = AtomIteratorFiltered.makeIterator(
                    new AtomIteratorSequenceDirectable(), typeFilter);
            ((AtomIteratorSequenceDirectable) inner.getWrappedIterator())
                    .setNumToSkip(1);
            return new ApiIntragroup(new ApiInnerVariable(outer,
                    (AtomIteratorAtomDependent) inner));
        }
        return makeIntergroupTypeIterator(types);
    }
}
