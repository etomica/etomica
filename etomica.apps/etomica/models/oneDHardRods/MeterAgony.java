package etomica.models.oneDHardRods;

import etomica.api.IPotential;
import etomica.data.DataSourceScalar;
import etomica.potential.PotentialMasterMonatomic;

public class MeterAgony extends DataSourceScalar {
    int numTrials, numAccept;
    IPotential potentialTarget, potentialHarmonic;
    PotentialMasterMonatomic pm;
    
    private static final long serialVersionUID = 1L;
    
    
    public MeterAgony(){
        
    }
    @Override
    public double getDataAsScalar() {
        // TODO Auto-generated method stub
        return 0;
    }

}
