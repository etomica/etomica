package etomica.modules.glass2d;

import etomica.atom.AtomType;
import etomica.atom.IAtomList;
import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.space.Vector;
import etomica.units.dimensions.CompoundDimension;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Length;
import etomica.units.dimensions.Time;

import java.util.Arrays;

public class DataSourceMSD implements IDataSource, ConfigurationStorage.ConfigurationStorageListener, DataSourceIndependent {

    protected final ConfigurationStorage configStorage;
    protected DataDoubleArray tData;
    protected DataDoubleArray.DataInfoDoubleArray tDataInfo;
    protected DataFunction data;
    protected DataFunction.DataInfoFunction dataInfo;
    protected double[] msdSum;
    protected final DataTag tTag, tag;
    protected long[] nSamples;
    protected final AtomType type;

    public DataSourceMSD(ConfigurationStorage configStorage) {
        this(configStorage, null);
    }

    public DataSourceMSD(ConfigurationStorage configStorage, AtomType type) {
        this.configStorage = configStorage;
        this.type = type;
        msdSum = new double[0];
        nSamples = new long[0];
        tag = new DataTag();
        tTag = new DataTag();
        reset();
    }

    public void reset() {
        int n = configStorage.getLastConfigIndex();
        if (n + 1 == msdSum.length && data != null) return;
        if (n < 1) n = 0;
        else n--;
        msdSum = Arrays.copyOf(msdSum, n);
        nSamples = Arrays.copyOf(nSamples, n);
        data = new DataFunction(new int[]{n});
        tData = new DataDoubleArray(new int[]{n});
        tDataInfo = new DataDoubleArray.DataInfoDoubleArray("t", Time.DIMENSION, new int[]{n});
        tDataInfo.addTag(tTag);
        dataInfo = new DataFunction.DataInfoFunction("MSD", new CompoundDimension(new Dimension[]{Length.DIMENSION}, new double[]{2}), this);
        dataInfo.addTag(tag);
        double[] t = tData.getData();
        if (t.length > 0) {
            double[] savedTimes = configStorage.getSavedTimes();
            double dt = savedTimes[0] - savedTimes[1];
            for (int i = 0; i < t.length; i++) {
                t[i] = dt * (1L << i);
            }
        }
    }

    @Override
    public IData getData() {
        if (configStorage.getLastConfigIndex() < 1) return data;
        double[] y = data.getData();
        for (int i = 0; i < msdSum.length; i++) {
            y[i] = msdSum[i] / nSamples[i];
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

    @Override
    public void newConfigruation() {
        reset(); // reallocates if needed
        long step = configStorage.getSavedSteps()[0];
        Vector[] positions = configStorage.getSavedConfig(0);
        IAtomList atoms = configStorage.getBox().getLeafList();
        for (int i = 1; i < msdSum.length; i++) {
            if (step % (1L << (i - 1)) == 0) {
                Vector[] iPositions = configStorage.getSavedConfig(i);
                for (int j = 0; j < positions.length; j++) {
                    if (type != null && atoms.get(j).getType() != type) continue;
                    msdSum[i - 1] += positions[j].Mv1Squared(iPositions[j]);
                    nSamples[i - 1]++;
                }
            }
        }
    }

    @Override
    public DataDoubleArray getIndependentData(int i) {
        return tData;
    }

    @Override
    public DataDoubleArray.DataInfoDoubleArray getIndependentDataInfo(int i) {
        return tDataInfo;
    }

    @Override
    public int getIndependentArrayDimension() {
        return 1;
    }

    @Override
    public DataTag getIndependentTag() {
        return tTag;
    }
}
