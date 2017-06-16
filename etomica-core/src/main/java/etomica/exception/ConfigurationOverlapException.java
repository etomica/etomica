/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.exception;

import etomica.box.Box;

/**
 * Exception thrown when an overlap is detected in the configuration of
 * a box, subject to the status of an ignoreOverlap flag.
 */
public class ConfigurationOverlapException extends RuntimeException {

    /**
     * Constructor for ConfigurationOverlapException
     * 
     * @param box the Box in which the overlap is detected
     */
    public ConfigurationOverlapException(Box box) {
        super("Overlap in configuration "+box);
        this.box = box;
    }

    private static final long serialVersionUID = 1L;
    public final Box box;
}
