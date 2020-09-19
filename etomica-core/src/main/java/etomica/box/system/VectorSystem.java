package etomica.box.system;

import etomica.space.Space;
import etomica.space.Vector;

import java.util.Arrays;

public class VectorSystem {
    private double[] data;
    private Vector[] viewVectors;
    private int count;
    private final int stride;
    private final Space space;


    public VectorSystem(Space space, int count) {
        this.stride = space.D();
        this.data = new double[count * stride];
        this.count = count;
        this.viewVectors = new Vector[count];
        this.space = space;
        this.makeVectors();
    }

    private void makeVectors() {
        for (int i = 0; i < viewVectors.length; i++) {
            switch (space.D()) {
                case 3:
                    viewVectors[i] = new ViewVector3D(i * 3, data);
                    break;
                default:
                    throw new IllegalArgumentException();
            }
        }
    }

    public Vector get(int i) {
        return this.viewVectors[i];
    }

    public void resize(int newSize) {
        this.count = newSize;
        this.data = Arrays.copyOf(this.data, newSize * stride);
        this.viewVectors = new Vector[newSize];
        this.makeVectors();
    }

    public int size() {
        return this.count;
    }
}
