/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import etomica.atom.AtomType;
import etomica.potential.IteratorDirective;

/**
 * Provides static methods for construction of some useful pair iterators.
 */
public final class ApiBuilder {

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
        AtomIteratorArrayListAdjacent aiInnerUp = new AtomIteratorArrayListAdjacent(IteratorDirective.Direction.UP);
        AtomIteratorArrayListAdjacent aiInnerDn = new AtomIteratorArrayListAdjacent(IteratorDirective.Direction.DOWN);
        return new ApiIntragroup(aiInnerUp, aiInnerDn);
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
        return makeNonAdjacentPairIterator(2);
    }

    public static ApiIntragroup makeNonAdjacentPairIterator(int numToSkip) {
        AtomIteratorArrayList aiInnerUp = new AtomIteratorArrayList(IteratorDirective.Direction.UP, numToSkip);
        AtomIteratorArrayList aiInnerDn = new AtomIteratorArrayList(IteratorDirective.Direction.DOWN, numToSkip);
        return new ApiIntragroup(aiInnerUp, aiInnerDn);
    }
        
        
    /**
     * Returns an intergroup iterator that filters the iterates so that only
     * those having the given type instances are returned. The given types are
     * applied corresponding to the pair of basis atoms identified when the
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
        if (types[0].getSpecies() == types[1].getSpecies() && types[0] != types[1]) {
            AtomIteratorBasisFilteredType2 outer = new AtomIteratorBasisFilteredType2(types[0], types[1]);
            AtomIteratorBasisFilteredType2 inner = new AtomIteratorBasisFilteredType2(types[0], types[1]);
            return new ApiIntergroupIntraSpecies(outer, inner);
        }
        AtomIteratorBasisFilteredType outer = new AtomIteratorBasisFilteredType(types[0]);
        AtomIteratorBasisFilteredType inner = new AtomIteratorBasisFilteredType(types[1]);
        return new ApiIntergroup(outer, inner);
    }
}
