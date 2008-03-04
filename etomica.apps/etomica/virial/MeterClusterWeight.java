package etomica.virial;

import etomica.api.IBox;
import etomica.api.IPotentialMaster;

import etomica.atom.iterator.IteratorDirective;
import etomica.data.DataSourceScalar;
import etomica.units.Dimension;
import etomica.units.Null;

/**
 * Meter to calculate the sampling weight of a cluster configuration.   
 */
 
public class MeterClusterWeight extends DataSourceScalar {
    
    public MeterClusterWeight(IPotentialMaster potentialMaster) {
        super("Cluster Weight",Null.DIMENSION);
        potential = potentialMaster;
    }
      
    public Dimension getDimension() {return Null.DIMENSION;}
    
    public double getDataAsScalar() {
    	weight.reset();
    	potential.calculate(box, up, weight);
    	return weight.weight();
    }

    /**
     * @return Returns the box.
     */
    public IBox getBox() {
        return box;
    }
    /**
     * @param box The box to set.
     */
    public void setBox(IBox box) {
        this.box = box;
    }

    private IBox box;
    private final IPotentialMaster potential;
    private final PotentialCalculationClusterWeightSum weight = new PotentialCalculationClusterWeightSum();
    private final IteratorDirective up = new IteratorDirective();
}