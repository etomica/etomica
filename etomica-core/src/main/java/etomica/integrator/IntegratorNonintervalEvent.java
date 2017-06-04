/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator;


/**
 * Event thrown by integrator when it announces reaching special points in the 
 * simulation process, such as its beginning and end. 
 *
 * @author David Kofke and Andrew Schultz
 *
 */
public class IntegratorNonintervalEvent extends IntegratorEvent {

    private static final long serialVersionUID = 1L;

    public IntegratorNonintervalEvent(Integrator source, NonintervalEventType type) {
        super(source, type);
    }

    public static final NonintervalEventType RESET = new NonintervalEventType("Reset"); //integrator is resetting

    public static class NonintervalEventType extends Type {
        private static final long serialVersionUID = 1L;

        protected NonintervalEventType(String label) {
            super(label);
        }
        
        public static Type[] choices() {
            return new NonintervalEventType[] {RESET};
        }
    }
}
