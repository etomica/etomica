/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator;

import etomica.util.EnumeratedType;

/**
 * Event object given to registered listeners when an integrator fires
 * any type event.  Object provides information about the type of event
 * being fired, and gives a reference to the integrator.
 *
 * @author David Kofke
 *
 */
public class IntegratorEvent {

    // Typed constants used to indicate the type of event integrator is
    // announcing
    
    protected final Type type;
    protected Integrator source;

    /**
     * 
     */
    public IntegratorEvent(Integrator source, Type type) {
        this.source = source;
        this.type = type;
    }

    public Type type() {
        return type;
    }
    
    public Integrator getSource() {
        return source;
    }

    //class used to mark the different types of interval events
    public static class Type extends EnumeratedType {
        private static final long serialVersionUID = 1L;

        protected Type(String label) {
            super(label);
        }

        public static Type[] choices() {
            return new Type[]{IntegratorNonintervalEvent.RESET};
        }
    }

}
