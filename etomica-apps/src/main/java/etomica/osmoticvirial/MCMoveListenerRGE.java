package etomica.osmoticvirial;

import etomica.atom.IAtom;
import etomica.atom.iterator.AtomIterator;
import etomica.box.Box;
import etomica.data.AccumulatorAverage;
import etomica.data.DataInfo;
import etomica.data.types.DataDoubleArray;
import etomica.integrator.mcmove.MCMove;
import etomica.integrator.mcmove.MCMoveEvent;
import etomica.integrator.mcmove.MCMoveTrialCompletedEvent;
import etomica.integrator.mcmove.MCMoveTrialFailedEvent;
import etomica.species.ISpecies;
import etomica.species.Species;
import etomica.units.dimensions.Null;
import etomica.util.IEvent;
import etomica.util.IListener;


public class MCMoveListenerRGE implements IListener {
    DataDoubleArray data;
    AccumulatorAverage accumulatorAverage;
    Box box;
    ISpecies species;
    int numAtoms;

    public MCMoveListenerRGE(AccumulatorAverage accumulatorAverage, Box box, ISpecies species, int numAtoms){
        DataInfo dataInfo = new DataDoubleArray.DataInfoDoubleArray("0/1", Null.DIMENSION, new int[]{numAtoms/2+1});
        accumulatorAverage.putDataInfo(dataInfo);
        this.accumulatorAverage = accumulatorAverage;
        this.box = box;
        this.species = species;
        this.numAtoms = numAtoms;
        this.data = new DataDoubleArray(numAtoms/2+1);
    }

    @Override
    public void actionPerformed(IEvent evt) {
        if (!(evt instanceof MCMoveTrialCompletedEvent)) {
            return;
        }

        double[] x = data.getData();
        int molecules = box.getNMolecules(species);
        int index;
        if(molecules < numAtoms/2+1) index = molecules;
        else index = numAtoms - molecules;
        for(int i=0; i<x.length; i++){
            x[i] = 0;
        }
        x[index] = 1;
        accumulatorAverage.putData(data);

    }
}
