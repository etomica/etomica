package etomica.modules.adsorption;

import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IVector;
import etomica.data.DataSourceScalar;
import etomica.units.Quantity;

public class MeterExcessAdsorbed extends DataSourceScalar {

    protected IBox box;
    protected int dim;
    protected double xMin, xMax;
    protected double pressure;
    protected final EOSSW eos;
    
    public MeterExcessAdsorbed(EOSSW eos) {
        super("excess adsorbed", Quantity.DIMENSION);
        this.eos = eos;
    }
    
    public void setRange(int dim, double xMin, double xMax) {
        this.dim = dim;
        this.xMin = xMin;
        this.xMax = xMax;
    }
    
    public void setBox(IBox box) {
        this.box = box;
    }
    
    public void setPressure(double pressure) {
        this.pressure = pressure;
    }

    public double getDataAsScalar() {
        IAtomList list = box.getLeafList();
        int n = 0;
        for (int i=0; i<list.getAtomCount(); i++) {
            IVector p = list.getAtom(i).getPosition();
            double x = p.getX(dim);
            if (x>xMin && x<xMax) n++;
        }
        double v = box.getBoundary().volume();
        double Ly = box.getBoundary().getBoxSize().getX(1);
        v *= (xMax-xMin)/Ly;
        double nHomogeneous = eos.rhoForPressure(pressure)*v;
        double nExPerArea = (n - nHomogeneous)/(v/Ly);
        return nExPerArea;
    }

}
