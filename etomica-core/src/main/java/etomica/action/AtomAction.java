/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.action;

import etomica.atom.IAtom;

/**
 * Interface for a class that can perform an action on an atom.
 */
public interface AtomAction {

    public void actionPerformed(IAtom atom);
}
