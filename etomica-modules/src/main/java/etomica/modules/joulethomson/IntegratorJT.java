/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.joulethomson;

import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorManagerMC;
import etomica.modifier.ModifierBoolean;
import etomica.util.random.IRandom;

public class IntegratorJT extends IntegratorManagerMC {
    
    IntegratorBox integratorNPH;
    IntegratorBox integratorNVE;
    int nveCount;
    boolean doNVE = true;
    boolean wasReset = false;
    
    public IntegratorJT(IRandom random, IntegratorBox nph, IntegratorBox nve) {
        super(random);
        integratorNPH = nph;
        integratorNVE = nve;
        addIntegrator(nph);
        addIntegrator(nve);
    }
    
    public void reset() {
        super.reset();
        nveCount = 100;
    }

    protected void doStepInternal() {
        if(random.nextDouble() < globalMoveProbability) {
            doGlobalMoves();
        } else {
            if (nveCount > 0) {
                integratorNVE.doStep();
                nveCount--;
                if (nveCount == 0) {
                    integratorNPH.reset();
                }
            }
            else {
                integratorNPH.doStep();
            }
        }
    }
    
    //inner class used to toggle between NVE and NPH ensembles
    public class EnsembleToggler implements ModifierBoolean {
        public void setBoolean(boolean doNVE) {
            if (doNVE) {
                if (nveCount == 0) {
                    nveCount = 100;
                }
            }
            else {
                nveCount = 0;
            }
        }
        public boolean getBoolean() {
            return nveCount > 0;
        }
    }
    
    
}
