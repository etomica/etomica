package etomica.box.storage;

import etomica.space.Space;
import etomica.space.Vector;

import java.util.Arrays;
import java.util.BitSet;

public class VectorStorage extends DoubleStructStorage<Vector> {
    private final Space space;

    public VectorStorage(Space space, int count) {
        super(space.D(), count, getVectorClass(space));
        this.space = space;
    }

    private static Class<? extends Vector> getVectorClass(Space space) {
        switch (space.D()) {
            case 3:
                return ViewVector3D.class;
            default:
                throw new IllegalStateException("Unexpected value: " + space.D());
        }
    }

    protected Vector makeView(int idx) {
        switch (space.D()) {
            case 3:
                return new ViewVector3D(idx * 3, data);
            default:
                throw new IllegalArgumentException();
        }
    }

    @Override
    protected void updateIndex(Vector view, int newIdx) {
        switch (space.D()) {
            case 3:
                ((ViewVector3D) view).setIndex(newIdx * stride);
                break;
            default:
                throw new IllegalArgumentException();
        }
    }

    @Override
    protected void updateData(Vector[] views, double[] newData) {
        for (Vector view : views) {
            if (view == null) { continue; }

            switch (space.D()) {
                case 3:
                    ((ViewVector3D) view).setData(newData);
                    break;
                default:
                    throw new IllegalArgumentException();
            }
        }
    }
}
