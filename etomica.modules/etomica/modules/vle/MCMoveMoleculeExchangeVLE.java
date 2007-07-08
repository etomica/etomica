package etomica.modules.vle;

import etomica.data.meter.MeterDensity;
import etomica.integrator.IntegratorBox;
import etomica.integrator.mcmove.MCMoveMoleculeExchange;
import etomica.potential.PotentialMaster;
import etomica.util.IRandom;

public class MCMoveMoleculeExchangeVLE extends MCMoveMoleculeExchange {

    public MCMoveMoleculeExchangeVLE(PotentialMaster potentialMaster, IRandom random, IntegratorBox integrator1, IntegratorBox integrator2) {
        super(potentialMaster, random, integrator1, integrator2);
        meterDensity = new MeterDensity(potentialMaster.getSpace());
    }
    
    public boolean doTrial() {
        boolean success = super.doTrial();
        if (!success) return success;
        meterDensity.setBox(box1);
        double density1 = meterDensity.getDataAsScalar();
        meterDensity.setBox(box2);
        double density2 = meterDensity.getDataAsScalar();
        if (density2 > density1) {
            rejectNotify();
            return false;
        }
        return true;
    }
    
    private static final long serialVersionUID = 1L;
    protected double density;
    protected MeterDensity meterDensity;
}
