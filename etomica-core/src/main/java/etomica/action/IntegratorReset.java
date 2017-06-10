/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.action;

import etomica.integrator.Integrator;
import etomica.exception.ConfigurationOverlapException;


/**
 * Action that calls the reset method of an integrator.
  */
public class IntegratorReset extends IntegratorActionAdapter {

    public IntegratorReset() {
        super();
    }
    
    public IntegratorReset(Integrator integrator, boolean ignoreOverlap) {
        super(integrator);
        this.ignoreOverlap = ignoreOverlap;
    }

    public void actionPerformed() {
        if(integrator != null) {
            try {
                integrator.reset();
            }
            catch (ConfigurationOverlapException e) {
                if (!ignoreOverlap) {
                    throw e;
                }
            }
        }
    }

    /**
     * @return Returns the ignoreOverlap.
     */
    public boolean isIgnoreOverlap() {
        return ignoreOverlap;
    }

    /**
     * @param ignoreOverlap The ignoreOverlap to set.
     */
    public void setIgnoreOverlap(boolean ignoreOverlap) {
        this.ignoreOverlap = ignoreOverlap;
    }

    private static final long serialVersionUID = 1L;
    private boolean ignoreOverlap;
}
