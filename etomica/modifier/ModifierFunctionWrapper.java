package etomica.modifier;

import etomica.util.Function;
import etomica.util.FunctionInvertible;

/**
 * Extends the configurable modifier to permit a function to be applied to
 * a value before the modified property is set.  Correspondingly, the value
 * obtained from a get of the modified property has applied the inverse
 * of the function.
 *
 * @see ModifierGeneral
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
     * Applies function to given value before setting modulated property.
     */
    public void setValue(double d) {
        super.setValue(function.f(d));
    }
    
    /**
     * Applies inverse of function to modulated property before returning it.
     */
    public double getValue() {
        return function.inverse(super.getValue());
    }
    
    private static final long serialVersionUID = 1L;
}