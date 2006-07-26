package etomica.models.hexane;

import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.DataTag;
import etomica.data.meter.Meter;
import etomica.data.types.DataGroup;
import etomica.data.types.DataTensor;
import etomica.data.types.DataGroup.DataInfoGroup;
import etomica.phase.Phase;

public class MeterNormalMode implements Meter {
       
    public MeterNormalMode(Phase ph, PairIndexerMolecule pri, MeterCorrelationMatrix mcm) {
        phase = ph;
        dim = phase.space().D();
        pi = pri;
        tag = new DataTag();
        
        //The data that is gotten is a series of flattened tensors that describes 
        //nan what does it describe?
        //he average of gamma
        dataCorMat = updateMCM(mcm);
        
        
        
        
        
        resetMeter();
    }
    
    
    
    
    public DataTag getTag() {
        return tag;
    }
    
    //nan fill this space!!
    public Data getData(){
        resetMeter();
        
    }
    
    //nan fill this space!!
    public void resetMeter(MeterCorrelationMatrix mcm){
        updateMCM(mcm);
        
    }
    
    public void updateMCM(MeterCorrelationMatrix mcm){
//      The data that is gotten is a series of flattened tensors that describes 
        //nan what does it describe?
        //he average of gamma
        dataCorMat = mcm.getData();
        
        //do something to reassemble the tensor here!
        
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
    
    
    
    //probably we will need a phase, a pair indexer, and a metercorrelationmatrix



    private DataGroup data;
    private DataInfoGroup dataInfo;
    private Data dataCorMat;
    private String name;
    private Phase phase;
    private final DataTag tag;
    private int dim;        //The dimension of the phase and space.
    private PairIndexerMolecule pi;     //Used to find the bin each thing is stored in.
}
