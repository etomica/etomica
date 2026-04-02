package etomica.data.histogram;

import etomica.math.DoubleRange;
import etomica.space3d.Vector3D;

import java.util.ArrayList;

public interface HistogramVector {

    public void addVector (Vector3D vector3D);

    public double[] getXHistogram();

    public ArrayList<Vector3D> getVectorHistogram();

    public long getCount();

    public void setNBins(int n);

    public void setRange(DoubleRange rangeX, DoubleRange rangeY, DoubleRange rangeZ);

    public DoubleRange[] getRangeVector();

    public void reset();
}
