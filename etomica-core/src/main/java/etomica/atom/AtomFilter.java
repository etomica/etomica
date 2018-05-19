/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import java.util.function.Predicate;


/**
 * Interface for a class that screens atoms according
 * to some criterion.
 */
public interface AtomFilter extends Predicate<IAtom> {

    /**
     * Returns true if atom is passes test of filter.
     */
    boolean test(IAtom a);

}
