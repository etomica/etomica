/**
 * 
 */
package etomica.data.meter;

import etomica.atom.AtomSet;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.IAtomKinetic;
import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.DataSource;
import etomica.data.DataTag;
import etomica.data.IDataInfo;
import etomica.data.types.DataVector;
import etomica.data.types.DataVector.DataInfoVector;
import etomica.box.Box;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.units.CompoundDimension;
import etomica.units.Dimension;
import etomica.units.Length;
import etomica.units.Mass;
import etomica.units.Time;

/**
 * Returns the instantaneous total center-of-mass momentum, summed over all
 * leaf atoms in a box.
 *
 */
public class MeterMomentumCOM implements DataSource, java.io.Serializable {

    public MeterMomentumCOM(Space space) {
        data = new DataVector(space);
        momentumSum = data.x;
        dataInfo = new DataInfoVector("COM momentum", new CompoundDimension(
                new Dimension[] {Mass.DIMENSION, Length.DIMENSION, Time.DIMENSION}, new double[] {1.,1.,-1.}),
                space);
        tag = new DataTag();
        dataInfo.addTag(tag);
        
    }

    /**
     * Returns the instantaneous total center-of-mass momentum over all atoms in the box.
     */
    public Data getData() {
        momentumSum.E(0.0);
        AtomSet leafList = box.getSpeciesMaster().getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
            double mass = ((AtomTypeLeaf)a.getType()).getMass();
            momentumSum.PEa1Tv1(mass,a.getVelocity());
        }
        return data;
    }
    
    /**
     * @return Returns the box.
     */
    public Box getBox() {
        return box;
    }
    /**
     * @param box The box to set.
     */
    public void setBox(Box box) {
        this.box = box;
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
    private Box box;
    private final IVector momentumSum;
    private final DataVector data;    
    private final DataInfo dataInfo;
    private String name;
    protected final DataTag tag;

}
