package etomica.nbr.list;

import etomica.box.Box;
import etomica.box.BoxAgentManager;

public class NeighborListAgentSource implements BoxAgentManager.BoxAgentSource<NeighborListManager> {
    protected PotentialMasterList potentialMaster;
    protected double range;

    public NeighborListAgentSource(double range) {
        this.range = range;
    }

    public void setRange(double newRange) {
        range = newRange;
    }

    public void setPotentialMaster(PotentialMasterList p){
        potentialMaster = p;
    }

    public NeighborListManager makeAgent(Box box) {
        return new NeighborListManager(potentialMaster, range, box);
    }

    public void releaseAgent(NeighborListManager object) {
        object.dispose();
    }
}
