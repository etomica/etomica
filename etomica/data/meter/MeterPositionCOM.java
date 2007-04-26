/**
 * 
 */
package etomica.data.meter;

import etomica.atom.AtomLeaf;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.DataSource;
import etomica.data.DataTag;
import etomica.data.IDataInfo;
import etomica.data.types.DataVector;
import etomica.data.types.DataVector.DataInfoVector;
import etomica.phase.Phase;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.units.Length;
import etomica.util.NameMaker;

/**
 * Returns the instantaneous center-of-mass position, summed over all
 * leaf atoms in a phase, dividing by the number of atoms.
 *
 */
public class MeterPositionCOM implements DataSource, java.io.Serializable {

    public MeterPositionCOM(Space space) {
        data = new DataVector(space);
        positionSum = data.x;
        dataInfo = new DataInfoVector("COM momentum", Length.DIMENSION, space);
        setName(NameMaker.makeName(this.getClass()));
        tag = new DataTag();
        dataInfo.addTag(tag);
        
    }

    /**
     * Returns the position of the center of mass of all atoms in the phase.
     */
    public Data getData() {
        iterator.reset();
        positionSum.E(0.0);
        double massSum = 0.0;
        for (AtomLeaf atom = (AtomLeaf)iterator.nextAtom(); atom != null;
             atom = (AtomLeaf)iterator.nextAtom()) {
            double mass = ((AtomTypeLeaf)atom.getType()).getMass();
            massSum += mass;
            positionSum.PEa1Tv1(mass,atom.getPosition());
        }
        positionSum.TE(1.0/massSum);
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
        iterator.setPhase(phase);
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
    private final AtomIteratorLeafAtoms iterator = new AtomIteratorLeafAtoms();
    private final IVector positionSum;
    private final DataVector data;    
    private final DataInfo dataInfo;
    private String name;
    protected final DataTag tag;

}
