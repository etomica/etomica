package etomica.experimental;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.space.Vector;
import etomica.space3d.Vector3D;
import jdk.incubator.vector.DoubleVector;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class VecSys3DAlt {
    private static final DoubleVector.DoubleSpecies S = DoubleVector.preferredSpecies();

    public final double[] xs;
    public final double[] ys;
    public final double[] zs;

    public VecSys3DAlt(int rows) {
        this.xs = new double[rows];
        this.ys = new double[rows];
        this.zs = new double[rows];

        if(rows % S.length() != 0) {
            throw new RuntimeException();
        }
    }

    public VecSys3DAlt(Box box) {
        this(box.getLeafList().size());
        IAtomList leafList = box.getLeafList();
        for (int i = 0; i < leafList.size(); i++) {
            IAtom atom = leafList.get(i);
            Vector pos = atom.getPosition();
            xs[i] = pos.getX(0);
            ys[i] = pos.getX(1);
            zs[i] = pos.getX(2);
        }
    }

    public List<Vector> toVectors() {
        List<Vector> list = new ArrayList<>(this.size());
        for (int i = 0; i < this.size(); i++) {
            list.add(new Vector3D(xs[i], ys[i], zs[i]));
        }
        return list;
    }

    public int size() {
        return xs.length;
    }

    public void setVector(int i, Vector v) {
        xs[i] = v.getX(0);
        ys[i] = v.getX(1);
        zs[i] = v.getX(2);
    }

    public void setAll(double x) {
        Arrays.fill(xs, x);
        Arrays.fill(ys, x);
        Arrays.fill(zs, x);
    }

    public void addSingle(int i, Vector x) {
        this.xs[i] += x.getX(0);
        this.ys[i] += x.getX(1);
        this.zs[i] += x.getX(2);
    }

    public void subSingle(int i, Vector x) {
        this.xs[i] -= x.getX(0);
        this.ys[i] -= x.getX(1);
        this.zs[i] -= x.getX(2);
    }

    public Vector3D diffSingle(int i, int j) {
        double dx = this.xs[i] - this.xs[j];
        double dy = this.ys[i] - this.ys[j];
        double dz = this.zs[i] - this.zs[j];
        return new Vector3D(dx, dy, dz);
    }

    public void addScaled(double scale, VecSys3DAlt other) {
        for (int i = 0; i < xs.length; i++) {
            this.xs[i] += other.xs[i] * scale;
            this.ys[i] += other.ys[i] * scale;
            this.zs[i] += other.zs[i] * scale;
        }
//        int i;
//        var scales = S.broadcast(scale);
//        for (i = 0; i < xs.length; i += S.length()) {
//            var x = S.fromArray(other.xs, i).fma(scales, S.fromArray(this.xs, i));
//            x.intoArray(this.xs, i);
//
//            var y = S.fromArray(other.ys, i).fma(scales, S.fromArray(this.ys, i));
//            y.intoArray(this.ys, i);
//
//            var z = S.fromArray(other.zs, i).fma(scales, S.fromArray(this.zs, i));
//            z.intoArray(this.zs, i);
//        }

        // tail

    }
}
