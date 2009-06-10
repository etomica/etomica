package etomica.exception;

import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.atom.Atom;

/**
 * Exception thrown when an overlap is detected in the configuration of
 * a box, subject to the status of an ignoreOverlap flag.
 */
public class ConfigurationOverlapException extends RuntimeException {

    /**
     * Constructor for ConfigurationOverlapException
     * 
     * @param box the Box in which the overlap is detected
     */
    public ConfigurationOverlapException(IBox box) {
        super("Overlap in configuration "+box);
        this.box = box;

        IAtomList list = box.getLeafList();
        for(int i = 0; i < 32; i++){
            System.out.println("atom " + i);
            for(int j = 0; j < 3; j++){
                System.out.println(       ((Atom)list.getAtom(i)).getPosition().x(j)     );
            }
        }
        throw new IllegalStateException("Dude!");
    }

    private static final long serialVersionUID = 1L;
    public final IBox box;
}
