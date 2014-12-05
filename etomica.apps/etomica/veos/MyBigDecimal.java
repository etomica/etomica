/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.veos;

import java.math.BigDecimal;
import java.math.MathContext;

/**
 * Data class to handle operations involving large numbers.  Encapsulates a BigDecimal,
 * and supports operations to add, multiply, etc. with other instances. 
 * 
 * @author kofke
 *
 */

public class MyBigDecimal implements IBigValue {

    private BigDecimal value;
    private static final MyBigDecimal ZERO = new MyBigDecimal(BigDecimal.ZERO);
    private static final MathContext outputPrecision = new MathContext(16);
    private final MathContext workingPrecision = new MathContext(0);
    
    /**
     * Default constructor sets value to 0.0
     */
    public MyBigDecimal() {
       this(ZERO);
    }
    
    public IBigValue getZero() {
        return ZERO;
    }

    /**
     * Constructs using the given value.
     */
    public MyBigDecimal(double value) {
        this.E(value);
    }
        
    /**
     * Copy constructor.
     */
    public MyBigDecimal(IBigValue v) {
        this.E(v);
    }
    
    public MyBigDecimal(BigDecimal v) {
        value = v;
    }
    
    public boolean isZero() {
        return value.signum() == 0;
    }
    
    /**
     * Returns the value held by the instance (in true form, not as log, and with sign applied).
     */
    public double value() {
        return value.doubleValue();
    }
    
    public boolean isPositive() {
        return value.signum() > 0;
    }
    
    /**
     * Returns logarithm of the value represented by the instance. Returns NaN if the
     * value is negative.
     */
    public double lnValue() {
        //return BigDecimalMath.log(value).doubleValue();
        BigDecimal temp = value.round(outputPrecision);
        int scale = temp.scale();
        double unscaledValue = temp.unscaledValue().doubleValue();
        return Math.log(unscaledValue) - scale * Math.log(10.);
    }
        
    /**
     * Equals (=) operation.  Assigns the given value to this instance, replacing current value.
     */
    public final void E(double value) {
        this.value = BigDecimal.valueOf(value);
    }
    
    /**
     * Equals (=) operation. Assigns given value (specified via another instance of LogValue) to this instance, 
     * replacing current value.
     */
    public final void E(IBigValue x) {
        this.value = ((MyBigDecimal)x).value;
    }
        
    /**
     * Plus-equals (+=) operation.  Replaces the current value with the one obtained by adding
     * the given value, considering signs of both.
     */
    public void PE(IBigValue x) {
        value = value.add(((MyBigDecimal)x).value,workingPrecision);
    }
        
    /**
     * Plus-equals (+=) a1 times v1.
     */
    public void PEa1Tv1(double a1, IBigValue iv1) {
        BigDecimal av = (((MyBigDecimal)iv1).value).multiply(BigDecimal.valueOf(a1),workingPrecision);
        value = value.add(av,workingPrecision);
    }
    
    /**
     * Times-equals (*=) operation.  Replaces the current value with the one obtained by multiplying
     * the given value, considering signs of both.
     */
    public void TE(IBigValue x) {
        value = value.multiply(((MyBigDecimal)x).value,workingPrecision);
    }
    
    /**
     * Divide-equals (/=) operation.  Replaces the current value with the one obtained by dividing
     * the given value, considering signs of both.
     */
    public void DE(IBigValue x) {
        value = value.divide(((MyBigDecimal)x).value,new MathContext(200));
    }
    
    /**
     * Times-equals (*=). Replaces current value by multiplying it by the given value.
     * @param a1
     */
    public void TE(double a1) {
        value = value.multiply(BigDecimal.valueOf(a1),workingPrecision);
    }
    
    private static void test(double a, double b) {
        IBigValue lna = new MyBigDecimal(a);
        MyBigDecimal lnb = new MyBigDecimal(b);
//        lna.PE(lnb);
//        lnb.PE(new LogValue(a));
//        System.out.println((a+b)+"  "+lna.value()+"  "+lnb.value());
        lna.TE(lnb);
        lnb.TE(new MyBigDecimal(a));
        System.out.println((a*b)+"  "+lna.value()+"  "+lnb.value());
    }
    
    public static void main(String[] arg) {
        double a = 1.34e-25;
        double b = 1/a;//3.0e+70;
        System.out.println(a+" "+Math.exp(Math.log(a)));
        System.out.println(b+" "+Math.exp(Math.log(b)));
        MyBigDecimal.test(+a,+b);
        MyBigDecimal.test(-a,+b);
        MyBigDecimal.test(+a,-b);
        MyBigDecimal.test(-a,-b);
        MyBigDecimal.test(+b,+a);
        MyBigDecimal.test(-b,+a);
        MyBigDecimal.test(+b,-a);
        MyBigDecimal.test(-b,-a);
        MyBigDecimal ad = new MyBigDecimal(a);
        System.out.println(Math.log(a)+", "+ad.lnValue());
        
        BigDecimal var1 = new BigDecimal(0.0);
        BigDecimal var2 = new BigDecimal(17.0);
        BigDecimal var3 = var1.add(var2);
        System.out.println(var3.doubleValue());
        
    }
}
