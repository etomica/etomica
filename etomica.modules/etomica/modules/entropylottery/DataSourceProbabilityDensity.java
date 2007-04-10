package etomica.modules.entropylottery;

import etomica.atom.AtomLeaf;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.Data;
import etomica.data.DataSource;
import etomica.data.DataTag;
import etomica.data.IDataInfo;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.integrator.IntegratorIntervalListener;
import etomica.integrator.IntegratorNonintervalEvent;
import etomica.integrator.IntegratorNonintervalListener;
import etomica.integrator.IntegratorPhase;
import etomica.phase.Phase;
import etomica.space.BoundaryPeriodic;
import etomica.space.IVector;
import etomica.units.Quantity;

public class DataSourceProbabilityDensity implements DataSource, IntegratorIntervalListener, IntegratorNonintervalListener {

    public DataSourceProbabilityDensity() {
        dataInfo = new DataInfoDoubleArray("probability density", Quantity.DIMENSION, new int[]{0});
        data = new DataDoubleArray(0);
        atomIterator = new AtomIteratorLeafAtoms();
        newData = new double[0];
        tag = new DataTag();
        dataInfo.addTag(tag);
    }
    
    public Data getData() {
        return data;
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    public int getPriority() {
        return 0;
    }

    public void intervalAction() {
        data.assignTo(newData);
        double[] oldData = data.getData();
        int nBin = oldData.length;
        if (((BoundaryPeriodic)phase.getBoundary()).getPeriodicity()[0]) {
            // with PBC, let balls drift around the boundary
            newData[0] += ((oldData[nBin-1]+oldData[1])/2 - oldData[0])/totalAtomCount;
            newData[nBin-1] += ((oldData[nBin-1]+oldData[0])/2 - oldData[nBin-1])/totalAtomCount;
        }
        else {
            newData[0] += (oldData[1]-oldData[0])/(2*totalAtomCount); 
            newData[nBin-1] += (oldData[nBin-2]-oldData[nBin-1])/(2*totalAtomCount);
        }
        for (int i=1; i<nBin-1; i++) {
            newData[i] += ((oldData[i-1]+oldData[i+1])/2 - oldData[i])/totalAtomCount;
        }
        data.E(newData);
    }

    public void nonintervalAction(IntegratorNonintervalEvent evt) {
        if (evt.type() == IntegratorNonintervalEvent.RESET) {
            phase = ((IntegratorPhase)evt.getSource()).getPhase();
            totalAtomCount = phase.getSpeciesMaster().moleculeCount();
            IVector dimensions = phase.getBoundary().getDimensions();
            if (data.getLength() != (int)Math.round(dimensions.x(0))) {
                int newSize = (int)Math.round(dimensions.x(0));
                data = new DataDoubleArray(newSize);
                dataInfo = new DataInfoDoubleArray("probability density", Quantity.DIMENSION, new int[]{newSize});
                dataInfo.addTag(tag);
                newData = new double[newSize];
            }
            data.E(0);
            atomIterator.setPhase(phase);
            atomIterator.reset();
            double[] atomCount = data.getData();
            while(atomIterator.hasNext()) {
                AtomLeaf a = (AtomLeaf)atomIterator.nextAtom();
                int x = (int)Math.round(a.getCoord().getPosition().x(0)+dimensions.x(0)*0.5-0.5);
                atomCount[x]++;
            }
        }
    }

    protected DataDoubleArray data;
    protected DataInfoDoubleArray dataInfo;
    protected Phase phase;
    protected final AtomIteratorLeafAtoms atomIterator;
    protected double[] newData;
    protected int totalAtomCount;
    protected final DataTag tag;
}
