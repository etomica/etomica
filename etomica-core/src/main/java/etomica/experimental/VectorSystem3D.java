package etomica.experimental;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.space.Vector;
import etomica.space2d.Vector2D;
import etomica.space3d.Vector3D;

import java.util.Arrays;

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
    
    public void setVector(int row, Vector v) {
        coords[row * 3] = v.getX(0);
        coords[row * 3 + 1] = v.getX(1);
        coords[row * 3 + 2] = v.getX(2);
    }

    public void setAll(double x) {
        Arrays.fill(this.coords, x);
    }

    public int getRows() {
        return this.rows;
    }

    public Vector3D get(int i) {
        return new Vector3D(coords[i * 3], coords[i * 3 + 1], coords[i * 3 + 2]);
    }
    
    public Vector3D diff(int v1, int v2) {
        double dx = coords[3 * v1] - coords[3 * v2];
        double dy = coords[3 * v1 + 1] - coords[3 * v2 + 1];
        double dz = coords[3 * v1 + 2] - coords[3 * v2 + 2];
        return new Vector3D(dx, dy, dz);
    }

    public void add(int i, Vector3D vec) {
        coords[3 * i] += vec.getX(0);
        coords[3 * i + 1] += vec.getX(1);
        coords[3 * i + 2] += vec.getX(2);
    }

    public void sub(int i, Vector3D vec) {
        coords[3 * i] -= vec.getX(0);
        coords[3 * i + 1] -= vec.getX(1);
        coords[3 * i + 2] -= vec.getX(2);
    }

    public void addScaled(int i, int j, double scale, VectorSystem3D sys) {
        this.coords[i * 3 + 0] += scale * sys.coords[j * 3 + 0];
        this.coords[i * 3 + 1] += scale * sys.coords[j * 3 + 1];
        this.coords[i * 3 + 2] += scale * sys.coords[j * 3 + 2];
    }
}
