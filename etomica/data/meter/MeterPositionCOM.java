/**
 * 
 */
package etomica.data.meter;

import etomica.atom.AtomSet;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.IAtomPositioned;
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
import etomica.units.Length;

/**
 * Returns the instantaneous center-of-mass position, summed over all
 * leaf atoms in a box, dividing by the number of atoms.
 *
 */
public class MeterPositionCOM implements DataSource, java.io.Serializable {

    public MeterPositionCOM(Space space) {
        data = new DataVector(space);
        positionSum = data.x;
        dataInfo = new DataInfoVector("COM momentum", Length.DIMENSION, space);
        tag = new DataTag();
        dataInfo.addTag(tag);
        
    }

    /**
     * Returns the position of the center of mass of all atoms in the box.
     */
    public Data getData() {
        positionSum.E(0.0);
        double massSum = 0.0;
        AtomSet leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomPositioned a = (IAtomPositioned)leafList.getAtom(iLeaf);
            double mass = ((AtomTypeLeaf)a.getType()).getMass();
            massSum += mass;
            positionSum.PEa1Tv1(mass,a.getPosition());
        }
        positionSum.TE(1.0/massSum);
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
    public void setBox(Box newBox) {
        box = newBox;
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
    private final IVector positionSum;
    private final DataVector data;    
    private final DataInfo dataInfo;
    private String name;
    protected final DataTag tag;

}
