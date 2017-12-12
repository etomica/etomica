package etomica.server.representations;

import etomica.space.Boundary;

public class ConfigurationUpdate {
    private final double[][][] coordinates;
    private final Boundary[] boxBoundaries;

    public ConfigurationUpdate(double[][][] coordinates, Boundary[] boxBoundaries) {
        this.coordinates = coordinates;
        this.boxBoundaries = boxBoundaries;
    }

    public Boundary[] getBoxBoundaries() {
        return boxBoundaries;
    }

    public double[][][] getCoordinates() {
        return coordinates;
    }
}
