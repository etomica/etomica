package etomica.osmoticvirial;

import etomica.box.Box;
import etomica.data.AccumulatorAverage;
import etomica.data.DataInfo;
import etomica.data.types.DataDoubleArray;
import etomica.species.Species;
import etomica.util.IEvent;
import etomica.util.IListener;


public class MCMoveListenerRGE implements IListener {
    DataDoubleArray data;
    AccumulatorAverage accumulatorAverage;
    Box box;
    Species species;
    int numAtoms;

    public MCMoveListenerRGE(AccumulatorAverage accumulatorAverage, Box box, Species species, int numAtoms){
        DataInfo dataInfo = new DataDoubleArray.DataInfoDoubleArray("0/1", null, new int[]{numAtoms+1});
        accumulatorAverage.putDataInfo(dataInfo);
        this.accumulatorAverage = accumulatorAverage;
        this.box = box;
        this.species = species;
        this.numAtoms = numAtoms;
        this.data = new DataDoubleArray(numAtoms+1);
    }

    @Override
    public void actionPerformed(IEvent event) {
        double[] x = data.getData();
        int index = box.getNMolecules(species);
        for(int i=0; i<x.length; i++){
            x[i] = 0;
        }
        x[index] = 1;
        accumulatorAverage.putData(data);
    }
}
