package etomica.models.oneDHardRods;

import etomica.api.IPotential;
import etomica.data.DataSourceScalar;
import etomica.potential.PotentialMasterMonatomic;
import etomica.units.Null;

public class MeterAgony extends DataSourceScalar {
    int numTrials, numAccept;
    IPotential potentialTarget, potentialHarmonic;
    PotentialMasterMonatomic pm;
    
    private static final long serialVersionUID = 1L;
    
    
    public MeterAgony(){
        super("agony", Null.DIMENSION);
    }
    @Override
    public double getDataAsScalar() {
        // TODO Auto-generated method stub
        return 0;
    }

}
