/**
 * 
 */
package etomica.data.meter;

import etomica.atom.AtomArrayList;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomTypeLeaf;
import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.DataSource;
import etomica.data.DataTag;
import etomica.data.IDataInfo;
import etomica.data.types.DataVector;
import etomica.data.types.DataVector.DataInfoVector;
import etomica.phase.Phase;
import etomica.space.ICoordinateKinetic;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.units.CompoundDimension;
import etomica.units.Dimension;
import etomica.units.Length;
import etomica.units.Mass;
import etomica.units.Time;
import etomica.util.NameMaker;

/**
 * Returns the instantaneous total center-of-mass momentum, summed over all
 * leaf atoms in a phase.
 *
 */
public class MeterMomentumCOM implements DataSource, java.io.Serializable {

    public MeterMomentumCOM(Space space) {
        data = new DataVector(space);
        momentumSum = data.x;
        dataInfo = new DataInfoVector("COM momentum", new CompoundDimension(
                new Dimension[] {Mass.DIMENSION, Length.DIMENSION, Time.DIMENSION}, new double[] {1.,1.,-1.}),
                space);
        setName(NameMaker.makeName(this.getClass()));
        tag = new DataTag();
        dataInfo.addTag(tag);
        
    }

    /**
     * Returns the instantaneous total center-of-mass momentum over all atoms in the phase.
     */
    public Data getData() {
        momentumSum.E(0.0);
        AtomArrayList leafList = phase.getSpeciesMaster().getLeafList();
        int nLeaf = leafList.size();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            AtomLeaf a = (AtomLeaf)leafList.get(iLeaf);
            double mass = ((AtomTypeLeaf)a.getType()).getMass();
            momentumSum.PEa1Tv1(mass,((ICoordinateKinetic)a).getVelocity());
        }
        return data;
    }
    
    /**
     * @return Returns the phase.
     */
    public Phase getPhase() {
        return phase;
    }
    /**
     * @param phase The phase to set.
     */
    public void setPhase(Phase phase) {
        this.phase = phase;
    }
    
    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }
    
    public DataTag getTag() {
        return tag;
    }
    
    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    private static final long serialVersionUID = 1L;
    private Phase phase;
    private final IVector momentumSum;
    private final DataVector data;    
    private final DataInfo dataInfo;
    private String name;
    protected final DataTag tag;

}
