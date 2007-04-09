package etomica.modules.joulethomson;

import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.IntegratorManagerMC;
import etomica.integrator.IntegratorPhase;
import etomica.modifier.ModifierBoolean;
import etomica.potential.PotentialMaster;
import etomica.util.IRandom;

public class IntegratorJT extends IntegratorManagerMC {
    
    IntegratorPhase integratorNPH;
    IntegratorPhase integratorNVE;
    int nveCount;
    boolean doNVE = true;
    boolean wasReset = false;
    
    public IntegratorJT(PotentialMaster potentialMaster, IRandom random, IntegratorPhase nph, IntegratorPhase nve) {
        super(potentialMaster, random);
        integratorNPH = nph;
        integratorNVE = nve;
        addIntegrator(nph);
        addIntegrator(nve);
    }
    
    public void reset() throws ConfigurationOverlapException {
        super.reset();
        nveCount = 100;
    }

    public void doStepInternal() {
        if(random.nextDouble() < globalMoveProbability) {
            doGlobalMoves();
        } else {
            if (nveCount > 0) {
                integratorNVE.doStep();
                nveCount--;
                if (nveCount == 0) {
                    try {
                        integratorNPH.reset();
                    }
                    catch (ConfigurationOverlapException e) {
                        // this shouldn't happen
                        throw new RuntimeException(e);
                    }
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