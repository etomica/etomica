package etomica;

/**
 * A consumer of data.  Data goes in and might (or might not) come out.
 */
public interface DataSink {
	public abstract void add(double[] values);
}