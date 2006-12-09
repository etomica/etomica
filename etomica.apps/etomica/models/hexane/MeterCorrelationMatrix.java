package etomica.models.hexane;

import etomica.atom.AtomLeaf;
import etomica.atom.AtomPair;
import etomica.atom.iterator.ApiLeafAtoms;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.DataTag;
import etomica.data.meter.Meter;
import etomica.data.types.DataGroup;
import etomica.data.types.DataTensor;
import etomica.data.types.DataGroup.DataInfoGroup;
import etomica.data.types.DataTensor.DataInfoTensor;
import etomica.phase.Phase;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.units.Area;
import etomica.units.Null;

/**
 * Calculates a correlation matrix.
 * @author nancycribbin
 *
 */
public class MeterCorrelationMatrix implements Meter {

    public MeterCorrelationMatrix(Phase ph, PairIndexerMolecule pi) {
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
        tag = new DataTag();
        dataInfo.addTag(tag);
        
        //Set up the pair iterator.
        api1 = new ApiLeafAtoms();
        api1.setPhase(phase);
        api1.reset();
        
        //Set up the single atom iterator.
        atomIterator.setPhase(phase);
        atomIterator.reset();
        op = new Vector[phase.getSpeciesMaster().getMaxGlobalIndex()+1];
        //set up the vectors where the original position is stored
        while (atomIterator.hasNext()) {
            AtomLeaf next = (AtomLeaf)atomIterator.nextAtom();
            op[next.getGlobalIndex()] = (Vector)next.getCoord().position().clone();
        }
        
        //Set up the tempVex.
        maxlength = pri.getMaxLength();
        int tempInt = pri.getMaxLength();
        tempVex = new Vector[tempInt];
        for(int i = 0; i < tempInt; i++){
            tempVex[i] = phase.space().makeVector();
        }
        counter = new int[maxlength];
        tempTen = phase.space().makeTensor();
        tempTenAgain = phase.space().makeTensor();
        
        resetMeter();
        setSymmetric(true);
    }
    
    public DataTag getTag() {
        return tag;
    }

    /**
     * returns the correlation matrix 
     */
    public Data getData() {
        resetMeter();
        
        //Used in a variety of loops; save some time by making them now
        AtomLeaf atom;
        int bin;
        Vector vex = phase.space().makeVector();
        
        if(!symmetric) System.out.println(symmetric);
        // make the change-in-position-from-the-original-lattice-point vector
        atomIterator.reset();
        while (atomIterator.hasNext()) {
            atom = (AtomLeaf) atomIterator.nextAtom();
            tempVex[atom.getGlobalIndex()].E(atom.getCoord().position());
            tempVex[atom.getGlobalIndex()].ME(op[atom.getGlobalIndex()]);
        }

        //Calculate the deltas tensor, and store it.
        api1.reset();        
        while (api1.hasNext()) {
            AtomPair pair = api1.nextPair();
            bin = pri.getBin(pair.getAtom(0), pair.getAtom(1));
            Vector gam1 = tempVex[pair.getAtom(0).getGlobalIndex()];
            Vector gam2 = tempVex[pair.getAtom(1).getGlobalIndex()];
            tempTen.Ev1v2(gam1, gam2);
            if(symmetric){
                tempTenAgain.E(tempTen);
                tempTen.transpose();
                tempTen.PE(tempTenAgain);
                ((DataTensor)data.getData(bin)).x.PE(tempTen);
                counter[bin] += 2;
            } else {
                ((DataTensor)data.getData(bin)).x.PE(tempTen);
                counter[bin] += 1;
            }
        }
        api1.reset();
        
        //Calculate the deltas tensor for each atom with itself.
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            atom = (AtomLeaf)atomIterator.nextAtom();
            bin = pri.getBin(atom, atom);
            vex.E(atom.getCoord().position());
            vex.ME(op[atom.getGlobalIndex()]);
            tempTen.Ev1v2(vex, vex);
//            System.out.println(atom.getGlobalIndex() +"  "+ bin +"  "+ vex.x(0)+"  "+ vex.x(1) +"  "+ vex.x(2));
            if(symmetric){
                tempTenAgain.E(tempTen);
                tempTen.transpose();
                tempTen.PE(tempTenAgain);
                ((DataTensor)data.getData(bin)).x.PE(tempTen);
                counter[bin] += 2;
            } else {
                ((DataTensor)data.getData(bin)).x.PE(tempTen);
                counter[bin] += 1;
            }
        }
            
            
            

        //Divide each bin by the amount of times we put something in i
        for(int i = 0; i < maxlength; i++){
            if(counter[i] != 0){
                (((DataTensor)data.getData(i)).x).TE(1/(double)counter[i]);
            }
        }
        
        return data;
    }

    /**
     * Resets the iterators, zeros the bin count, and zeros each tensor.
     * 
     */
    public void resetMeter() {

        // Reset the array of tensors to zero
        api1.reset();
        atomIterator.reset();
        AtomLeaf firstAtom = (AtomLeaf) atomIterator.nextAtom();
        int ticker = 0;
        while (atomIterator.hasNext()) {
            int bin = pri.getBin(firstAtom, atomIterator.nextAtom());
            ticker++;
            ((DataTensor)data.getData(bin)).x.E(0.0);
        }
        for(int i = 0; i < maxlength; i++){
            counter[i] = 0;
        }
        tempTen.E(0.0);
        tempTenAgain.E(0.0);
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

    public boolean isSymmetric() {
        return symmetric;
    }

    public void setSymmetric(boolean symmetric) {
        this.symmetric = symmetric;
    }

    private PairIndexerMolecule pri; // Calculates the index to store stuff under
    private DataGroup data;
    private DataInfoGroup dataInfo;
    private String name;
    private Phase phase;
    int dim; // The dimension of the space.
    // This is the massive iterator which iterates over all the leaf atoms of
    // the phase. It returns pairs of leaf atoms.
    private final ApiLeafAtoms api1;
    private final AtomIteratorLeafAtoms atomIterator = new AtomIteratorLeafAtoms();
    Vector[] op; // The original positions of the atoms.
    private Vector[] tempVex;
    private final DataTag tag;
    private int maxlength;  //The maximum length of various things
    private int[] counter;
    private boolean symmetric;
    Tensor tempTen;
    Tensor tempTenAgain;

}
