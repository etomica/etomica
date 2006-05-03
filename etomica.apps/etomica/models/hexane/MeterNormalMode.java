package etomica.models.hexane;

import etomica.atom.AtomLeaf;
import etomica.atom.AtomPair;
import etomica.atom.iterator.ApiLeafAtoms;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.meter.Meter;
import etomica.data.types.DataGroup;
import etomica.data.types.DataTensor;
import etomica.data.types.DataGroup.DataInfoGroup;
import etomica.data.types.DataTensor.DataInfoTensor;
import etomica.phase.Phase;
import etomica.space.Vector;
import etomica.units.Area;
import etomica.units.Null;

public class MeterNormalMode implements Meter {

    public MeterNormalMode(Phase ph, PairIndexer pi) {
        this.phase = ph;
        dim = phase.space().D();
        this.pri = pi;

        DataTensor[] dataTensors = new DataTensor[pri.getMaxLength()];
        DataInfoTensor[] dataInfoTensors = new  DataInfoTensor[pri.getMaxLength()];
        for (int i = 0; i < pri.getMaxLength(); i++) {
            dataTensors[i] = new DataTensor(phase.space());
            dataInfoTensors[i] = new DataInfoTensor("Normal Mode deltas", Area.DIMENSION, phase.space());
        }
        data = new DataGroup(dataTensors);
        dataInfo = new DataInfoGroup("Normal Mode deltas", Null.DIMENSION, dataInfoTensors);
        
        //Set up the pair iterator.
        api1 = new ApiLeafAtoms();
        api1.setPhase(phase);
        api1.reset();
        
        //set up the vectors where the original position is stored
        for (int i = 0; i < pri.getMaxLength(); i++) {
            op[i] = phase.space().makeVector();
        }
        while (api1.hasNext()) {
            op[tempAtom.getGlobalIndex()].E(tempAtom.coord.position());
        }
        
        count = new int[pri.getMaxLength()];

        tempAtom = new AtomLeaf(phase.space());
        //tempAtom = (AtomLeaf) api1.peek();

        resetMeter();

    }

    public Data getData() {
        resetMeter();

        // make the change-in-position-from-the-original-lattice-point vector
        leafIter.reset();
        while (leafIter.hasNext()) {
            AtomLeaf atom = (AtomLeaf) leafIter.nextAtom();
            tempVex[atom.getGlobalIndex()].E(atom.coord.position());
            tempVex[atom.getGlobalIndex()].ME(op[atom.getGlobalIndex()]);
        }

        //Calculate the deltas tensor, and store it.
        api1.reset();        
        while (api1.hasNext()) {
            AtomPair pair = api1.nextPair();
            
            // Now we do the funky multiplication-to-get-the-tensor thing.
            for (int i = 0; i < pri.getMaxLength(); i++) {
                for (int j = 0; j < pri.getMaxLength(); j++) {
                    int bin = pri.getBin(pair);
                    ((DataTensor)data.getData(bin)).x.PEv1v2(tempVex[i], tempVex[j]);
                    count[bin] += 1;
                }
            }
        }
        api1.reset();

        // Divide each tensor by the number of times something is stored in it.

        for (int i = 0; i < pri.getMaxLength(); i++) {
            if (count[i] != 0) {
                //Divide the tensor by count here.
                ((DataTensor)data.getData(i)).x.TE(1.0/(double)count[i]);
            }
        }
        return data;
    }

    /**
     * Resets the iterators, zeros the bin count, and zeros each tensor.
     * 
     */
    public void resetMeter() {
        // Reset the bin counters to zero
        for (int i = 0; i < pri.getMaxLength(); i++) {
            count[i] = 0;
        }

        // Reset the array of tensors to zero
        api1.reset();
        tempAtom = (AtomLeaf) leafIter.peek();
        while (leafIter.hasNext()) {
            ((DataTensor)data.getData(tempAtom.getGlobalIndex())).x.E(0.0);
        }
        
        // Reset the iterators
        leafIter.reset();
        api1.reset();
    }

    public DataInfo getDataInfo() {
        return dataInfo;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getName() {
        return name;
    }

    public void setPhase(Phase phase) {
        this.phase = phase;
    }

    public Phase getPhase() {
        return phase;
    }


    private PairIndexer pri; // Calculates the index to store stuff under
    private DataGroup data;
    private DataInfoGroup dataInfo;
    private String name;
    private int[] count; // Stores the number of times a storage location is used.
    private Phase phase;
    int dim; // The dimension of the space.
    // This is the massive iterator which iterates over all the leaf atoms of
    // the phase. It returns pairs of leaf atoms.
    private final ApiLeafAtoms api1;
    private final AtomIteratorLeafAtoms leafIter = new AtomIteratorLeafAtoms();
    Vector[] op; // The original positions of the atoms.
    private AtomLeaf tempAtom;
    private Vector[] tempVex;
}
