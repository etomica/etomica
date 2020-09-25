package etomica.box.storage;

import etomica.space.IOrientation;
import etomica.space.Space;
import etomica.space3d.Orientation3D;
import etomica.space3d.OrientationFull3D;

import java.util.Arrays;

public class OrientationStorage {
    private double[] data;
    private IOrientation[] viewOrientations;
    private int count;
    private final int stride;
    private final Space space;
    private final boolean isAxisSymmetric;


    public OrientationStorage(Space space, int count, boolean isAxisSymmetric) {
        int stride = space.D();
        if (isAxisSymmetric) {
            stride *= 2;
        }
        this.stride = stride;
        this.data = new double[count * stride];
        this.count = count;
        this.viewOrientations = new IOrientation[count];
        this.space = space;
        this.isAxisSymmetric = isAxisSymmetric;
        this.makeVectors();
    }

    private void makeVectors() {
        for (int i = 0; i < viewOrientations.length; i++) {
            switch (space.D()) {
                case 3:
                    if (this.isAxisSymmetric) {
                        viewOrientations[i] = new OrientationFull3D(space,
                                new ViewVector3D(i * stride, this.data),
                                new ViewVector3D(i * stride + 3, this.data));
                    } else {
                        viewOrientations[i] = new Orientation3D(new ViewVector3D(i * stride, this.data));
                    }
                    break;
                default:
                    throw new IllegalArgumentException();
            }
        }
    }

    public IOrientation get(int i) {
        return this.viewOrientations[i];
    }

    public void resize(int newSize) {
        this.count = newSize;
        this.data = Arrays.copyOf(this.data, newSize * stride);
        this.viewOrientations = new IOrientation[newSize];
        this.makeVectors();
    }

    public int size() {
        return this.count;
    }
}
