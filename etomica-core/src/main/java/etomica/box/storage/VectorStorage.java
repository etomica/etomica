package etomica.box.storage;

import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Vector3D;

public class VectorStorage extends DoubleStructStorage<Vector> {
    private final Space space;
    private final int D;

    public VectorStorage(Space space) {
        this(space, 0);
    }

    public VectorStorage(Space space, int count) {
        super(space.D(), getVectorClass(space));
        this.addNull(count);
        this.space = space;
        this.D = space.D();
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

    public Space getSpace() {
        return this.space;
    }

    public Vector load(int i) {
        Vector v = space.makeVector();
        switch (space.D()) {
            case 3:
                v.setX(0, this.data[i * stride]);
                v.setX(1, this.data[i * stride + 1]);
                v.setX(2, this.data[i * stride + 2]);
                break;
            default:
                throw new IllegalArgumentException();
        }
        return v;
    }

    public static Vector diff(VectorStorage vectorStorage, int idx1, int idx2) {
        switch (vectorStorage.D) {
            case 3:
                Vector3D v = new Vector3D(
                        vectorStorage.data[idx1 * vectorStorage.stride] - vectorStorage.data[idx2 * vectorStorage.stride],
                        vectorStorage.data[idx1 * vectorStorage.stride + 1] - vectorStorage.data[idx2 * vectorStorage.stride + 1],
                        vectorStorage.data[idx1 * vectorStorage.stride + 2] - vectorStorage.data[idx2 * vectorStorage.stride + 2]
                );
                return v;
            default:
                throw new IllegalArgumentException();
        }
    }
}
