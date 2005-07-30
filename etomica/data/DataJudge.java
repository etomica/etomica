package etomica.data;

import etomica.Data;
import etomica.data.types.DataArithmetic;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataInteger;


/**
 * Assigns a ranking to a Data object to expess the extent that it meets criteria
 * specified by a DataSink or some other client.  A lower value indicates a better
 * fit, with zero the smallest assigned value (perfect fit).
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on Jul 24, 2005 by kofke
 */
public interface DataJudge {

    public int rankData(Data data, int index, int depth);

    /**
     * Returns true if extraction is directed to recurse into sub-DataGroups; returns false otherwise.
     */
    public boolean doRecurse();

   
    ////////// INNER CLASSES /////////////

    public static class Group implements DataJudge {
        
        public Group(DataJudge[] judges) {
            this.judges = (DataJudge[])judges.clone();
        }
        public int rankData(Data data, int index, int depth) {
            int rank = 0;
            for(int i=0; i<judges.length; i++) {
                rank += judges[i].rankData(data, index, depth);
            }
            return rank;
        }
        /**
         * Returns true.
         */
        public boolean doRecurse() {
            return true;
        }
    
        private final DataJudge[] judges;
    }

    public static class ByIndex implements DataJudge {
        final int index;
        public ByIndex(int index) {
            this.index = index;
        }
        /**
         * If rank matches, returns depth, otherwise returns MAX_RANK.  By returning
         * depth, preference is given to top-level Data instances at the requested rank.
         */
        public int rankData(Data data, int index, int depth) {
            if(this.index == index) return depth;
            return MAX_RANK;
        }
        /**
         * Returns false.
         */
        public boolean doRecurse() {
            return false;
        }
    }

    public static class ByDepth implements DataJudge {
        final int depth;
        public ByDepth(int depth) {
            this.depth= depth;
        }
        public int rankData(Data data, int index, int depth) {
            if(this.depth== depth) return 0;
            return MAX_RANK;
        }
        /**
         * Returns true.
         */
        public boolean doRecurse() {
            return true;
        }
    }
    
    public static class ByClass implements DataJudge {
        final Class targetClass;
        final boolean doRecurse;
        public ByClass(Class targetClass, boolean doRecurse) {
            if(targetClass == null) throw new NullPointerException("Must specify a non-null class");
            this.targetClass = targetClass;
            this.doRecurse = doRecurse;
        }
        public int rankData(Data data, int index, int depth) {
            if(targetClass == data.getClass()) return 0;
            return MAX_RANK;
        }
        public boolean doRecurse() {
            return doRecurse;
        }

    }

    public static final DataJudge INSTANCEOF_ARITHMETIC = new DataJudge() {
        public int rankData(Data data, int index, int depth) {
            if(data instanceof DataArithmetic) {
                return 10*depth;
            }
            return MAX_RANK;
        };
        
        /**
         * Returns true.
         */
        public boolean doRecurse() {
            return true;
        }
    };
    
    public static final DataJudge ASSIGNABLE_TO_DOUBLE = new DataJudge() {
        public int rankData(Data data, int index, int depth) {
            if(data instanceof DataDouble) {
                return 10*depth;
            } else if(data instanceof DataInteger) {
                return 1 + 10*depth;
            } else if(data instanceof DataDoubleArray && ((DataDoubleArray)data).getLength() == 1) {
                return 2 + 10*depth;
            } else if(data instanceof DataFunction && ((DataFunction)data).getLength() == 1) {
                return 3 + 10*depth;
            } else {
                return MAX_RANK;
            }
        }
        public boolean doRecurse() {
            return true;
        }
    };

    public static final int MAX_RANK = 1000;
    
}

