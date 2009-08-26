package etomica.modules.mu;

import etomica.api.IAtomPositioned;
import etomica.api.IBox;
import etomica.api.IMoleculeList;
import etomica.api.ISpecies;
import etomica.api.IVector;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.units.CompoundDimension;
import etomica.units.Dimension;
import etomica.units.Quantity;
import etomica.units.Volume;

public class MeterDensitySides implements IEtomicaDataSource {

    public MeterDensitySides(IBox box, ISpecies species) {
        this.box = box;
        this.species = species;
        data = new DataDoubleArray(2);
        tag = new DataTag();
        dataInfo = new DataInfoDoubleArray("density", new CompoundDimension(new Dimension[]{Quantity.DIMENSION, Volume.DIMENSION}, new double[]{1,-1}), new int[]{2});
        dataInfo.addTag(tag);
    }
    
    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    public IData getData() {
        data.E(0);
        double[] x = data.getData();
        IMoleculeList moleculeList = box.getMoleculeList(species);
        for (int i=0; i<moleculeList.getMoleculeCount(); i++) {
            IVector p = ((IAtomPositioned)moleculeList.getMolecule(i).getChildList().getAtom(0)).getPosition();
            if (p.getX(0) < 0) {
                x[0]++;
            }
            else {
                x[1]++;
            }
        }
        data.TE(2.0/box.getBoundary().volume());
        return data;
    }

    protected final IBox box;
    protected final ISpecies species;
    protected final DataTag tag;
    protected final DataDoubleArray data;
    protected final DataInfoDoubleArray dataInfo;
}
