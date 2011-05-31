package etomica.virial;

import etomica.api.IBox;
import etomica.api.IMoleculeList;
import etomica.api.IPotentialMolecular;

public class MayerFunctionMolecularThreeBody extends MayerFunctionThreeBody {

    protected final IPotentialMolecular p3;

    public MayerFunctionMolecularThreeBody(IPotentialMolecular p3) {
        this.p3 = p3;
    }

    protected double energy(IMoleculeList molecules, double[] r2) {
        return p3.energy(molecules);
    }

    public void setBox(IBox box) {
        p3.setBox(box);
        super.setBox(box);
    }
}
