package etomica.virial;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleBDArray;
import etomica.units.dimensions.Null;

import java.math.BigDecimal;
import java.math.MathContext;

public class MeterPercolation implements IDataSource {

    private final DataDoubleArray data;
    private final IDataInfo dataInfo;
    private final DataTag tag;
    private final Box box;
    private final double initialCoordinates[];

    /**
     * Constructor for MeterVirial.
     */
    public MeterPercolation(Box box) {
        this.box = box;
        initialCoordinates = new double[getBox().getLeafList().size()*getBox().getSpace().D()];
        for (int i=0; i<getBox().getLeafList().size(); i++){
            for (int j=0; j<getBox().getSpace().D(); j++){
                initialCoordinates[i * getBox().getSpace().D() + j] = getBox().getLeafList().get(i).getPosition().getX(j);
            }
        }
        data = new DataDoubleArray(getBox().getLeafList().size()*getBox().getSpace().D());
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("Cluster Value", Null.DIMENSION, new int[]{box.getLeafList().size()*box.getSpace().D()});
        tag = new DataTag();
        dataInfo.addTag(tag);
    }

    @Override
    public IData getData() {

        double x[] = data.getData();
        IAtomList list = box.getLeafList();
        for (int i=0; i<list.size(); i++){
            for (int j=0; j<box.getSpace().D(); j++){
                double dx = list.get(i).getPosition().getX(j) - initialCoordinates[i * getBox().getSpace().D() + j];
                x[i * getBox().getSpace().D() + j] = dx * dx;
            }
        }
        return data;
    }

    @Override
    public DataTag getTag() {
        return tag;
    }

    @Override
    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public Box getBox() {
        return box;
    }
}
