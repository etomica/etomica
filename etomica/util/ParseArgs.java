/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.util;

import java.io.File;
import java.lang.reflect.Field;

/**
 * Class that handles parsing of commandline arguments.  Arguments are used to
 * assign fields in the parameterWrapper, much as is done in ReadParameters for
 * reading input files.
 * 
 * @author Andrew Schultz
 */
public class ParseArgs {

    public ParseArgs() {
    }
    
    public ParseArgs(ParameterBase parameterWrapper) {
        setParameterWrapper(parameterWrapper);
    }
    
    public static void doParseArgs(ParameterBase parameterWrapper, String[] args) {
        ParseArgs parser = new ParseArgs(parameterWrapper);
        parser.parseArgs(args, true);
    }
    
    /**
     * Returns the parameter wrapper.
     */
    public ParameterBase getParameterWrapper() {
        return wrapper;
    }

    /**
     * Sets the parameterWrapper
     */
    public void setParameterWrapper(ParameterBase newParameterWrapper) {
        wrapper = newParameterWrapper;
        fields = wrapper.getClass().getFields();
    }

    /**
     * Parses each argument and attempts to match it with a field from
     * parameterWrapper and sets the field to the value from the next argument.
     * If the next argument is another option (or if there are no more
     * options), the value is taken to be 'true' with the hope that the field
     * is a boolean.
     * Options are accepted in the form "-numSteps" or "--numSteps".
     * This routine handles boolean, int, long, double, String.  Arrays are
     * handled if they are quoted.
     * 
     * -slanty -temperature 1 -alpha "0.3 0.5"
     *  ==> slanty = true;
     *  ==> temperature = 1;
     *  ==> alpha = {0.3, 0.5};
     */
    public void parseArgs(String[] args) {
        parseArgs(args, false);
    }
    
    public void parseArgs(String[] args, boolean firstArgFile) {
        if (args.length == 0) return;
        if (firstArgFile) {
            if (new File(args[0]).exists()) {
                ReadParameters paramReader = new ReadParameters(args[0], wrapper);
                paramReader.readParameters();
                args = (String[])Arrays.removeObject(args, args[0]);
            }
        }
        if (args.length == 1 && (args[0].equals("-help") || args[0].equals("-h"))) {
            String strOut = "options and defaults: ";
            int len = strOut.length();
            System.out.print(strOut);
            for (int i=0; i<fields.length; i++) {
                try {
                    strOut = " -"+fields[i].getName()+" "+fields[i].get(wrapper).toString();
                }
                catch (IllegalAccessException e) {}
                len += strOut.length();
                if (len > 80) {
                    len = strOut.length()+5;
                    strOut = "\n     "+strOut;
                }
                System.out.print(strOut);
            }
            System.out.println();
            System.exit(0);
        }
        for (int i=0; i<args.length; i++) {
            if (args[i].charAt(0) != '-') {
                throw new RuntimeException("encountered "+args[i]+" when I was expecting an option");
            }
            String token = args[i].replaceAll("^--?", "");
            String value = "";
            if (i+1 == args.length || (args[i+1].charAt(0) == '-' && !Character.isDigit(args[i+1].charAt(1)))) {
                // last argument or next argument is another option.  hope this is a boolean parameter
                value = "true";
            }
            else {
                value = args[i+1];
                i++;
            }

            boolean foundField = false;
            for (int j=0; j<fields.length; j++) {
                if (token.equals(fields[j].getName())) {
                    wrapper.setValue(fields[j],value);
                    foundField = true;
                    break;
                }
                else if (token.equalsIgnoreCase("no"+fields[j].getName()) && value.equals("true")) {
                    wrapper.setValue(fields[j], "false");
                    foundField = true;
                    break;
                }
            }
            if (!foundField) {
                throw new RuntimeException("don't know what to do with token: "+token);
            }
        }
    }
    
    protected ParameterBase wrapper;
    protected Field[] fields;
}
