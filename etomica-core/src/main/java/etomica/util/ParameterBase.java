/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.util;

import java.lang.reflect.Array;
import java.lang.reflect.Field;


/**
 * Base class for all parameter wrapper implementations.  Subclasses should
 * simply contain fields corresponding to each parameter.  If any parameters
 * are more complex than can be handled by this class, ParamParsers should be
 * added to this object.
 */
public abstract class ParameterBase {
    
    public ParameterBase() {
        parsers = new ParameterParser[0];
    }
    
    public void addParser(ParameterParser newParser) {
        parsers = (ParameterParser[])Arrays.addObject(parsers, newParser);
    }
    
    public void removeParser(ParameterParser oldParser) {
        parsers = (ParameterParser[])Arrays.removeObject(parsers, oldParser);
    }
    
    public ParameterParser[] getParsers() {
        return parsers;
    }
    
    public final void setValue(Field field, String value) {
        for (int i=0; i<parsers.length; i++) {
            try {
                parsers[i].setValue(this, field, value);
                return;
            }
            catch (IllegalArgumentException e) {
                // that parser doesn't handle that field, move on
            }
        }
        Class type = field.getType();
        try {
            if (type == String.class) {
                field.set(this,value);
            }
            else if (type == int.class) {
                field.setInt(this,Integer.parseInt(value));
            }
            else if (type == long.class) {
                field.setLong(this,Long.parseLong(value));
            }
            else if (type == float.class) {
                field.setFloat(this,Float.parseFloat(value));
            }
            else if (type == double.class) {
                field.setDouble(this,Double.parseDouble(value));
            }
            else if (type == boolean.class) {
                field.setBoolean(this, Boolean.parseBoolean(value));
            }
            else if (type == char.class) {
                field.setChar(this,value.charAt(0));
            }
            else if (type.isEnum()) {
                field.set(this, Enum.valueOf(type, value));
            }
            else if (type.isArray()) {
                Class subType = type.getComponentType();
                String[] strings = value.split(" +");
                if (strings.length == 1 && strings[0].length() == 0) {
                    // split gives us a single empty string instead of
                    // an empty array
                    strings = new String[0];
                }
                if (subType == int.class) {
                    int[] array = new int[strings.length];
                    for (int i=0; i<array.length; i++) {
                        array[i] = Integer.parseInt(strings[i]);
                    }
                    field.set(this,array);
                }
                else if (subType == long.class) {
                    long[] array = new long[strings.length];
                    for (int i=0; i<array.length; i++) {
                        array[i] = Long.parseLong(strings[i]);
                    }
                    field.set(this,array);
                }
                else if (subType == float.class) {
                    float[] array = new float[strings.length];
                    for (int i=0; i<array.length; i++) {
                        array[i] = Float.parseFloat(strings[i]);
                    }
                    field.set(this,array);
                }
                else if (subType == double.class) {
                    double[] array = new double[strings.length];
                    for (int i=0; i<array.length; i++) {
                        array[i] = Double.parseDouble(strings[i]);
                    }
                    field.set(this,array);
                }
                else if (subType == boolean.class) {
                    boolean[] array = new boolean[strings.length];
                    for (int i=0; i<array.length; i++) {
                        array[i] = Boolean.parseBoolean(strings[i]);
                    }
                    field.set(this,array);
                }
                else if (subType == char.class) {
                    char[] array = new char[strings.length];
                    for (int i=0; i<array.length; i++) {
                        array[i] = strings[i].charAt(0);
                    }
                    field.set(this,array);
                }
                else if (subType.isEnum()) {
                    Enum[] array = (Enum[])Array.newInstance(subType, strings.length);
                    for (int i=0; i<array.length; i++) {
                        array[i] = Enum.valueOf(subType, strings[i]);
                    }
                    field.set(this,array);
                }
                else if (subType == String.class) {
                    field.set(this,strings);
                }
                else {
                    throw new RuntimeException("Unrecognized type "+type.getName());
                }
            }
            else {
                throw new RuntimeException("Unrecognized type "+type.getName());
            }
        }
        catch (NumberFormatException ex) {
            // catches exceptions thrown by Type.parseType for
            // int, long, float and double
            throw new RuntimeException("error parsing "+value+" for "+field.getName(), ex);
        }
        catch (IllegalAccessException e) {
            throw new RuntimeException("Illegal access exception trying to set "+field.getName()+" to "+value, e);
        }
    }

    /**
     * Actually writes field name and value to a file using fileWriter.
     * If the type of field is unrecognized, the wrapper object itself will
     * be asked to write the field.
     */
    public final String getValueString(Field field) {
        for (int i=0; i<parsers.length; i++) {
            try {
                String value = parsers[i].getValueString(this, field);
                return value;
            }
            catch (IllegalArgumentException e) {
                // that parser doesn't handle that field, move on
            }
        }
        Class type = field.getType();
        try {
            if (type.isArray()) {
                Class subType = type.getComponentType();
                String value = "";
                if (subType == int.class) {
                    int[] array = (int[])field.get(this);
                    for (int i=0; i<array.length; i++) {
                        if (i==0) value = Integer.toString(array[i]);
                        else value += " " + array[i];
                    }
                }
                if (subType == long.class) {
                    long[] array = (long[])field.get(this);
                    for (int i=0; i<array.length; i++) {
                        if (i == 0) value = Long.toString(array[i]);
                        else value += " " + array[i];
                    }
                }
                else if (subType == double.class) {
                    double[] array = (double[])field.get(this);
                    for (int i=0; i<array.length; i++) {
                        if (i == 0) value = Double.toString(array[i]);
                        else value += " " + array[i];
                    }
                }
                else if (subType == boolean.class) {
                    boolean[] array = (boolean[])field.get(this);
                    for (int i=0; i<array.length; i++) {
                        if (i == 0) value = Boolean.toString(array[i]);
                        else value += " " + array[i];
                    }
                }
                else {
                    value = field.get(this).toString();
                }
                return value;
            }
            return field.get(this).toString();
        }
        catch (IllegalAccessException e) {
            throw new RuntimeException("Illegal access exception trying to get "+field.getName()+" from "+this, e);
        }
    }
    
    protected ParameterParser[] parsers;
}
