package etomica;

/**
 * Class for constructing linked lists of Potentials.
 * Each Linker points to one potential and another Linker, the next one in the list.
 *
 * @author David Kofke
 */
public class PotentialLinker implements java.io.Serializable {
    public final Potential potential;
    public PotentialLinker next;
    //Constructors
//    public PotentialLinker(Potential a) {potential = a;}
    public PotentialLinker(Potential a, PotentialLinker l) {potential = a; next = l;}
}
