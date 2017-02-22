/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.util;

import java.lang.reflect.Field;

public interface ParameterParser {

    /**
     * Sets obj's field to appropriate value based on the given string, if this
     * parser is appropriate for the given field, or throws an exception if
     * it's not.
     */
    public void setValue(Object obj, Field field, String value) throws IllegalArgumentException;

    /**
     * Returns the value of obj's field as a string if this parser is 
     * appropirate for the given field, or throws an exception if it's not.
     */
    public String getValueString(Object obj, Field field) throws IllegalArgumentException;

}