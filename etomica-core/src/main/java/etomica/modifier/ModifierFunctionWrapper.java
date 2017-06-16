/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modifier;

import etomica.math.function.Function;
import etomica.math.function.FunctionInvertible;

/**
 * Extends the configurable modifier to permit a function to be applied to
 * a value before the modified property is set.  Correspondingly, the value
 * obtained from a get of the property has applied the inverse of the function.
 *
 * @author David Kofke
 */
 
public class ModifierFunctionWrapper extends ModifierGeneral {
    
    private FunctionInvertible function = new Function.Identity();
    
    public ModifierFunctionWrapper(Object[] obj, String prop) {
        super(obj, prop);
    }
    public ModifierFunctionWrapper(Object obj, String prop) {
        super(obj, prop);
    }

    public void setFunction(FunctionInvertible f) {function = f;}
    public FunctionInvertible getFunction() {return function;}
    
    /**
     * Applies function to given value before setting modified property.
     */
    public void setValue(double d) {
        super.setValue(function.f(d));
    }
    
    /**
     * Applies inverse of function to property before returning it.
     */
    public double getValue() {
        return function.inverse(super.getValue());
    }
    
    private static final long serialVersionUID = 1L;
}
