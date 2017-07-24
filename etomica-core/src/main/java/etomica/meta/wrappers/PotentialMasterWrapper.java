package etomica.meta.wrappers;

import etomica.meta.SimulationModel;
import etomica.potential.PotentialMaster;

public class PotentialMasterWrapper extends Wrapper<PotentialMaster> {


    public PotentialMasterWrapper(PotentialMaster wrapped, SimulationModel simModel) {
        super(wrapped, simModel);
    }

}
