package etomica.modules.interfacial;

import etomica.api.IAtom;
import etomica.api.IAtomPositioned;
import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IData;
import etomica.api.IMolecule;
import etomica.api.INearestImageTransformer;
import etomica.api.IVector;
import etomica.data.DataSourceAtomic;
import etomica.data.DataTag;
import etomica.data.IEtomicaDataInfo;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.space.ISpace;
import etomica.units.Angle;

/**
 * Meter for collecting the molecular orientation of the dimer.  The value
 * returned is cos(theta), where theta is the angle the dimer makes with the
 * x axis.
 */
public class MeterOrientation implements DataSourceAtomic {
    
    public MeterOrientation(ISpace space) {
        dataInfo = new DataInfoDouble("orientation", Angle.DIMENSION);
        data = new DataDouble();
        tag = new DataTag();
        dr = space.makeVector();
    }
    
    public DataTag getTag() {
        return tag;
    }
    
    public void setBox(IBox newBox) {
        transformer = newBox.getBoundary();
    }
    
    public IData getData(IAtom atom) {
        IAtomList children = ((IMolecule)atom).getChildList();
        dr.Ev1Mv2(((IAtomPositioned)children.getAtom(children.getAtomCount()-1)).getPosition(),
                  ((IAtomPositioned)children.getAtom(0)).getPosition());
        transformer.nearestImage(dr);
        data.x= dr.x(0) / Math.sqrt(dr.squared());
        return data;
    }
    
    public IEtomicaDataInfo getAtomDataInfo() {
        return dataInfo;
    }
    
    private static final long serialVersionUID = 1L;
    protected final DataInfoDouble dataInfo;
    protected final DataDouble data;
    protected final DataTag tag;
    protected final IVector dr;
    protected INearestImageTransformer transformer;
}
