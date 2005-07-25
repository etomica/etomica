package etomica.virial;

/**
 * @author andrew
 *
 * Clusters implementing ClusterWeight must return positive "values"
 */
public interface ClusterWeight extends ClusterAbstract {

    public interface Factory {
        public ClusterWeight makeWeightCluster(ClusterAbstract[] clusters);
    }
    
}
