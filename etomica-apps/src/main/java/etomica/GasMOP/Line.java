package etomica.GasMOP;

import etomica.space.Vector;
import etomica.space3d.Vector3D;

public class Line {
    public Vector point;
    public Vector direction;

    public Line(Vector point, Vector direction) {
        this.point = point;
        this.direction = direction;
    }
}
