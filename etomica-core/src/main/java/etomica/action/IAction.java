/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.action;

/**
 * General interface for performing an action, allowing the defined behavior to be separated from the driver that causes
 * it.  Actions may, for example, be invoked as listeners to events, or in response to user input.
 * Alternatively, they may be defined as a convenient encapsulation of nontrivial behaviors that are performed in
 * different contexts.
 */
public interface IAction {

    /**
     * Completes the action defined by the class implementing this interface.
     */
    void actionPerformed();

}