/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.util;

import etomica.api.IFunction;

/**
 * Interface for the basic features of a functional, which maps a function onto a double.
 */

public interface Functional {
    
    public double f(IFunction x);
 
}