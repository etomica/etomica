package etomica.data.types;


/**
 * Interface for a class that is capable of being converted into a DataArithmetic instance.
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on Jun 17, 2005 by kofke
 */
public interface DataNumeric {
    public DataArithmetic toArithmetic(DataArithmetic data);
}
