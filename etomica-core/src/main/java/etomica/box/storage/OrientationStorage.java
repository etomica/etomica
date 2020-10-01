package etomica.box.storage;

import etomica.space.IOrientation;
import etomica.space.Space;
import etomica.space2d.Orientation2D;
import etomica.space3d.Orientation3D;
import etomica.space3d.OrientationFull3D;

public class OrientationStorage extends DoubleStructStorage<IOrientation> {
    private final Space space;
    private final boolean isAxisSymmetric;


    public OrientationStorage(Space space, int count, boolean isAxisSymmetric) {
        super(isAxisSymmetric ? space.D() : space.D() * 2, count, IOrientation.class);
        this.space = space;
        this.isAxisSymmetric = isAxisSymmetric;
    }

    protected IOrientation makeView(int i) {
        switch (space.D()) {
            case 1:
                throw new IllegalArgumentException("We don't know how to make 1D Orientation that uses storage");
            case 2:
                return new Orientation2D(new ViewVector2D(i * stride, this.data));
            case 3:
                if (!this.isAxisSymmetric) {
                    return new OrientationFull3D(space,
                            new ViewVector3D(i * stride, this.data),
                            new ViewVector3D(i * stride + 3, this.data));
                } else {
                    return new Orientation3D(new ViewVector3D(i * stride, this.data));
                }
            default:
                throw new IllegalArgumentException();
        }
    }

    @Override
    protected void updateIndex(IOrientation view, int newIdx) {
        switch (space.D()) {
            case 1:
                ((ViewVector1D) view.getDirection()).setIndex(newIdx * stride);
                break;
            case 2:
                ((ViewVector2D) view.getDirection()).setIndex(newIdx * stride);
                break;
            case 3:
                ((ViewVector3D) view.getDirection()).setIndex(newIdx * stride);
                if (!this.isAxisSymmetric) {
                    OrientationFull3D o = (OrientationFull3D) view;
                    ((ViewVector3D) o.getSecondaryDirection()).setIndex(newIdx * stride + 3);
                }
                break;
            default:
                throw new IllegalStateException();
        }
    }

    @Override
    protected void updateData(IOrientation[] views, double[] newData) {
        for (IOrientation view : views) {
            if (view == null) {
                continue;
            }

            switch (space.D()) {
                case 1:
                    ((ViewVector1D) view.getDirection()).setData(newData);
                    break;
                case 2:
                    ((ViewVector2D) view.getDirection()).setData(newData);
                    break;
                case 3:
                    ((ViewVector3D) view.getDirection()).setData(newData);
                    if (!this.isAxisSymmetric) {
                        OrientationFull3D o = (OrientationFull3D) view;
                        ((ViewVector3D) o.getSecondaryDirection()).setData(newData);
                    }
                    break;
                default:
                    throw new IllegalStateException();
            }
        }
    }
}
