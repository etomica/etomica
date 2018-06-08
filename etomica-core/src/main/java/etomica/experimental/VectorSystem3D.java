package etomica.experimental;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.space.Vector;

public class VectorSystem3D {
    private final double[] coords;
    private final int rows;

    public VectorSystem3D(int vectors) {
        this.coords = new double[vectors * 3];
        this.rows = vectors;
    }

    public VectorSystem3D(Box box) {
        this(box.getLeafList().size());
        IAtomList leafList = box.getLeafList();
        for (int i = 0; i < leafList.size(); i++) {
            IAtom atom = leafList.get(i);
            Vector pos = atom.getPosition();
            coords[i * 3] = pos.getX(0);
            coords[i * 3 + 1] = pos.getX(1);
            coords[i * 3 + 2] = pos.getX(2);
        }
    }
}
