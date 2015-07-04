package etomica.normalmode;
import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IVectorMutable;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.space.ISpace;
import etomica.units.Null;

/**
 * Meter uses harmonically-mapped averages to compute properties and also
 * couples configurations with their mirror images to obtain cancellation
 * of error at very low temperatures. 
 */
 
public class MeterSolidMirror implements IEtomicaDataSource {
    
    protected IEtomicaDataSource ds;
    protected final DataTag tag;
    protected DataInfoDoubleArray dataInfo;
    protected DataDoubleArray data;
    protected double temperature;
    protected final IBox box;
    protected final CoordinateDefinition coordinateDefinition;
    protected final IVectorMutable dr;

    public MeterSolidMirror(ISpace space, IEtomicaDataSource dataSource, CoordinateDefinition coordinateDefinition) {
        this.coordinateDefinition = coordinateDefinition;
        this.ds = dataSource;
        tag = new DataTag();
        IEtomicaDataInfo di = dataSource.getDataInfo();
        box = coordinateDefinition.getBox();
        
        int n = di.getLength();
        dataInfo = new DataInfoDoubleArray("Stuff", Null.DIMENSION, new int[]{n});
        dataInfo.addTag(tag);
        data = new DataDoubleArray(n);
        dr = space.makeVector();
    }
    
    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }
    
    public DataTag getTag() {
        return tag;
    }
    
    public void setTemperature(double temperature) {
        this.temperature = temperature;
    }

    protected void reflect() {
        IAtomList atoms = box.getLeafList();
        for (int i=0; i<atoms.getAtomCount(); i++) {
            IAtom atom = atoms.getAtom(i);
            dr.Ev1Mv2(atom.getPosition(), coordinateDefinition.getLatticePosition(atom));
            atom.getPosition().PEa1Tv1(-2, dr);
        }
    }
    
    /**
     * Computes total pressure in box by summing virial over all pairs, and adding
     * ideal-gas contribution.
     */
    public IData getData() {
        data.E(ds.getData());
        double u0 = data.getValue(4);
        reflect();
        IData d1 = ds.getData();
        reflect();
        double u1 = d1.getValue(4);
        double f1 = Math.exp(-(u1-u0)/temperature);
        if (f1 < 1e-20) {
            return data;
        }
        if (f1 > 1e20) {
            data.E(d1);
            return data;
        }
        data.TE(1/(1+f1));
        d1.TE(f1/(1+f1));
        data.PE(d1);
        return data;
    }

}
