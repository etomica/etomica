package etomica;

import etomica.utility.Function;

/**
 * Extends the configurable modulator to permit a function to be applied to
 * a value before the modulated property is set.  Correspondingly, the value
 * obtained from a get of the modulated property has applied the inverse
 * of the function.
 *
 * @see Modulator
 * @author David Kofke
 */
 
public class ModulatorFunction extends Modulator {
    
    private Function function = new Function.Identity();
    
    public ModulatorFunction(Object[] obj, String prop) {
        super(obj, prop);
    }
    public ModulatorFunction(Object obj, String prop) {
        super(obj, prop);
    }

    public void setFunction(Function f) {function = f;}
    public Function getFunction() {return function;}
    
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
}