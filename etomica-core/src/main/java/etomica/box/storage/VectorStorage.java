package etomica.box.storage;

import etomica.space.Space;
import etomica.space.Vector;

public class VectorStorage extends DoubleStructStorage<Vector> {
    private final Space space;

    public VectorStorage(Space space) {
        this(space, 0);
    }

    public VectorStorage(Space space, int count) {
        super(space.D(), getVectorClass(space));
        this.addNull(count);
        this.space = space;
    }

    private static Class<? extends Vector> getVectorClass(Space space) {
        switch (space.D()) {
            case 1:
                return ViewVector1D.class;
            case 2:
                return ViewVector2D.class;
            case 3:
                return ViewVector3D.class;
            default:
                throw new IllegalStateException("Unexpected value: " + space.D());
        }
    }

    protected Vector createObject(int idx) {
        switch (space.D()) {
            case 1:
                return new ViewVector1D(idx, data);
            case 2:
                return new ViewVector2D(idx * 2, data);
            case 3:
                return new ViewVector3D(idx * 3, data);
            default:
                throw new IllegalArgumentException();
        }
    }

    @Override
    protected void updateIndex(Vector view, int newIdx) {
        switch (space.D()) {
            case 1:
                ((ViewVector1D) view).setIndex(newIdx * stride);
                break;
            case 2:
                ((ViewVector2D) view).setIndex(newIdx * stride);
                break;
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
                case 1:
                    ((ViewVector1D) view).setData(newData);
                    break;
                case 2:
                    ((ViewVector2D) view).setData(newData);
                    break;
                case 3:
                    ((ViewVector3D) view).setData(newData);
                    break;
                default:
                    throw new IllegalArgumentException();
            }
        }
    }
}
