package etomica.modules.droplet;

import etomica.api.IAtomList;
import etomica.api.IAtomPositioned;
import etomica.api.IBox;
import etomica.api.IData;
import etomica.api.IVectorMutable;
import etomica.data.DataTag;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.space.ISpace;
import etomica.space.Tensor;
import etomica.units.Null;

public class MeterDeformation implements IEtomicaDataSource {

    public MeterDeformation(ISpace space) {
        data = new DataDoubleArray(2);
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("deformation", Null.DIMENSION, new int[]{2});
        tag = new DataTag();
        f = new double[10];
        for (int i=0; i<f.length; i++) {
            f[i] = 1 + i/3.0;
        }
        
        center = space.makeVector();
        dr = space.makeVector();
        moment = space.makeTensor();
        workTensor = space.makeTensor();
    }
    
    public void setBox(IBox newBox) {
        box = newBox;
    }
    
    public IBox getBox() {
        return box;
    }
    
    public IData getData() {
        IAtomList leafList = box.getLeafList();
        for (int i=0; i<leafList.getAtomCount(); i++) {
            center.PE(((IAtomPositioned)leafList.getAtom(i)).getPosition());
        }
        center.TE(1.0/leafList.getAtomCount());
        
        for (int i=0; i<leafList.getAtomCount(); i++) {
            dr.Ev1Mv2(((IAtomPositioned)leafList.getAtom(i)).getPosition(), center);
            workTensor.Ev1v2(dr, dr);
            moment.PE(workTensor);
        }
        moment.TE(1.0/leafList.getAtomCount());
        
        double b1 = moment.trace();
        
        double b2 = 0;
        workTensor.E(moment);
        workTensor.TE(moment);
        for (int i=0; i<moment.D(); i++) {
            for (int j=0; j<moment.D(); j++) {
                b2 += workTensor.component(i,j);
            }
        }
        
        double b3 = moment.determinant();
        
        double c1 = b1 / (3.0*Math.pow(b3, 1.0/3.0)) - 1.0;
        double c2 = b2 / (3.0*Math.pow(b3, 2.0/3.0)) - 1.0;

        double delta = c2 - 4.0*c1;
        double factor = Math.signum(delta);

        double eOld = 0;
        double eOldOld = 0;
        double e = factor * Math.sqrt(c1/3.0);
        while (Math.abs(e-eOld) > 1.e-9 && (Math.abs(e-eOld) <= Math.abs(e-eOldOld) || (e-eOld)*(e-eOldOld) < 0)) {
            eOldOld = eOld;
            eOld = e;
            double denominator = 1.0;
            for (int i=0; i<10; i++) {
                denominator += f[i] * Math.pow(eOld, i);
            }
            e = factor * Math.sqrt(c1 / denominator);
        }
        
        double[] x = data.getData();
        x[0] = Math.pow(125*b3, 1.0/6.0);
        if (e > 0.02) {
            double a = Math.pow(1-e, 1.5);
            x[1] = (1.0 - a) / (1 + a);
        }
        else {
            double e2 = e*e;
            x[1] = 0.75*e + 0.375*e2 + 7.0/64.0*e2*e - 3.0/128.0*e2*e2 + 33.0/512.0*e2*e2*e - 61.0/1024.0*e2*e2*e2; 
        }
        return data;
    }

    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        // TODO Auto-generated method stub
        return null;
    }

    protected IBox box;
    protected final DataDoubleArray data;
    protected final DataDoubleArray.DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
    protected final double[] f;
    protected final IVectorMutable center;
    protected final IVectorMutable dr;
    protected final Tensor moment, workTensor;
}
