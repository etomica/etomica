/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.veos;

/**
 * Data class to handle operations involving large numbers.  Stores a value as the logarithm of its magnitude,
 * and supports operations to add, multiply, etc. with other instances. Handles positive and negative values
 * by keeping a separate field indicating sign of value.  
 * 
 * References to "value" in the comments indicate the value in its true form, not as logarithm.
 * 
 * @author kofke
 *
 */

public class BigValue implements IBigValue {

    private double lnValue;
    private boolean isPositive;
    public static final BigValue ZERO = new BigValue(0.0);
    
    /**
     * Default constructor sets value to 0.0
     */
    public BigValue() {
       this(0.0);
    }

    /**
     * Constructs using the given value.
     */
    public BigValue(double value) {
        this.E(value);
    }
    
    /**
     * Constructs given ln(|value|) and sign of value.
     */
    public BigValue(double lnValue, boolean isPos) {
        this.lnValue = lnValue;
        isPositive = isPos;
    }
    
    /**
     * Copy constructor.
     */
    public BigValue(IBigValue v) {
        this.E(v);
    }
    
    public IBigValue getZero() {
        return ZERO;
    }
    
    public boolean isZero() {
        return lnValue == Double.NEGATIVE_INFINITY;
    }
    
    /**
     * Returns the value held by the instance (in true form, not as log, and with sign applied).
     */
    public double value() {
        if(this.isZero()) return 0.0;
        return (isPositive ? +1 : -1) * Math.exp(lnValue);
    }
    
    public boolean isPositive() {
        return isPositive;
    }
    
    /**
     * Returns logarithm of the value represented by the instance. Returns NaN if the
     * value is negative.
     */
    public double lnValue() {
        if(isPositive) return lnValue;
        return Double.NaN;
    }
    
    /**
     * Equals (=) operation.  Assigns the given value to this instance, replacing current value.
     */
    public final void E(double value) {
        if(value == 0.0) {
            lnValue = Double.NEGATIVE_INFINITY;
            isPositive = false;
        } else {
            lnValue = Math.log(Math.abs(value));
            isPositive = (value > 0);
        }
    }
    
    /**
     * Equals (=) operation. Assigns given value (specified via another instance of LogValue) to this instance, 
     * replacing current value.
     */
    public final void E(IBigValue lv) {
        this.lnValue = ((BigValue)lv).lnValue;
        this.isPositive = ((BigValue)lv).isPositive;
    }
        
    /**
     * Plus-equals (+=) operation.  Replaces the current value with the one obtained by adding
     * the given value, considering signs of both.
     */
    public void PE(IBigValue ilnx) {
        BigValue lnx = (BigValue)ilnx;
        if(lnx.isZero()) return;
        if(this.isZero()) {
            this.E(lnx);
        } else {
            this.PE(lnx.lnValue, lnx.isPositive);
        }
    }
    
    private void PE(double lnxlnValue, boolean lnxisPositive) {
        //both have same sign
        if(this.isPositive == lnxisPositive) {
            if(this.lnValue > lnxlnValue) {
                lnValue += Math.log(1 + Math.exp(lnxlnValue-lnValue));
            } else {
                lnValue = lnxlnValue + Math.log(1 + Math.exp(lnValue-lnxlnValue)); 
            }
        } else { //different signs
            if(this.lnValue > lnxlnValue) {
                lnValue += Math.log(1 - Math.exp(lnxlnValue-lnValue));
            } else {
                lnValue = lnxlnValue + Math.log(1 - Math.exp(lnValue-lnxlnValue));
                isPositive = lnxisPositive;
            }
        }
    }
    
    /**
     * Plus-equals (+=) a1 times v1.
     */
    public void PEa1Tv1(double a1, IBigValue iv1) {
        BigValue v1 = (BigValue)iv1;
        if(a1 == 0 || v1.isZero()) return;
        if(this.isZero()) {
            this.E(v1);
            this.TE(a1);
        }
        boolean isPositiveSave = isPositive;
        double lnValueSave = lnValue;
        this.E(v1);
        this.TE(a1);
        this.PE(lnValueSave, isPositiveSave);
    }
    
    /**
     * Times-equals (*=) operation.  Replaces the current value with the one obtained by multiplying
     * the given value, considering signs of both.
     */
    public void TE(IBigValue ix) {
        BigValue x = (BigValue)ix;
        if(x.isZero()) {
            this.E(0.0);
        } else if(!this.isZero()) {
            lnValue += x.lnValue;
            if(!x.isPositive) isPositive = !isPositive;
        }
    }
    
    /**
     * Divide-equals (/=) operation.  Replaces the current value with the one obtained by dividing
     * the given value, considering signs of both.
     */
    public void DE(IBigValue ix) {
        BigValue x = (BigValue)ix;
        if(x.isZero()) {
            lnValue = Double.NaN;
            isPositive = false;
        } else if(!this.isZero()) {
            lnValue -= x.lnValue;
            if(!x.isPositive) isPositive = !isPositive;
        }
    }
    
    /**
     * Times-equals (*=). Replaces current value by multiplying it by the given value.
     * @param a1
     */
    public void TE(double a1) {
        if(a1 == 0) {
            this.E(0.0);
        } else if(!this.isZero()) {
            if(a1 > 0) {
                this.TEln(Math.log(a1));
            } else if(a1 < 0) {
                this.TEln(Math.log(-a1));
                isPositive = !isPositive;
            }
        }
    }

    
    /**
     * Times-equals (*=) operation, multiplying by a value given as its logarithm.  Value (i.e., the value
     * whose log is given as the argument) is assumed to be positive.
     */
    public void TEln(double lnx) {
        lnValue += lnx;
    }
    
    private static void test(double a, double b) {
        IBigValue lna = new BigValue(a);
        BigValue lnb = new BigValue(b);
//        lna.PE(lnb);
//        lnb.PE(new LogValue(a));
//        System.out.println((a+b)+"  "+lna.value()+"  "+lnb.value());
        lna.TE(lnb);
        lnb.TE(new BigValue(a));
        System.out.println((a*b)+"  "+lna.value()+"  "+lnb.value());
    }
    
    public static void main(String[] arg) {
        double a = 7.4e-15;
        double b = 1/a;//3.0e+70;
        System.out.println(a+" "+Math.exp(Math.log(a)));
        System.out.println(b+" "+Math.exp(Math.log(b)));
        BigValue.test(+a,+b);
        BigValue.test(-a,+b);
        BigValue.test(+a,-b);
        BigValue.test(-a,-b);
        BigValue.test(+b,+a);
        BigValue.test(-b,+a);
        BigValue.test(+b,-a);
        BigValue.test(-b,-a);
    }
}
