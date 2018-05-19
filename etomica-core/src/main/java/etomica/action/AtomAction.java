/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.action;

import etomica.atom.IAtom;

/**
 * Performs an action on an atom. Does not extend IAction because the actionPerformed method takes an IAtom argument.
 */
public interface AtomAction {

    /**
     * Performs the defined action on the specified atom
     * @param atom the specified atom
     */
    public void actionPerformed(IAtom atom);
}
