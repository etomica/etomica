package etomica.modules.entropylottery;

import etomica.action.AtomActionTranslateTo;
import etomica.api.IVector;
import etomica.atom.IAtom;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.config.Configuration;
import etomica.box.Box;
import etomica.space.Space;

/**
 * Configuration that simply puts all the Atoms at 0, or -0.5 if the box size
 * is even (which makes the entropy lottery binning happy).
 * @author Andrew Schultz
 */
public class ConfigurationZero implements Configuration, java.io.Serializable {

	private final Space space;

    public ConfigurationZero(Space _space) {
        super();
        this.space = _space;
    }

    public void initializeCoordinates(Box box) {
        AtomActionTranslateTo atomActionTranslateTo = new AtomActionTranslateTo(space);
        IVector work = space.makeVector();
        work.E(0.0);
        int intD = (int)Math.round(box.getBoundary().getDimensions().x(0));
        if (intD % 2 == 0) {
            work.E(-0.5);
        }
        atomActionTranslateTo.setDestination(work);

        AtomIteratorAllMolecules iterator = new AtomIteratorAllMolecules(box);
        iterator.reset();
        for (IAtom a = iterator.nextAtom(); a != null; a = iterator.nextAtom()) {
           atomActionTranslateTo.actionPerformed(iterator.nextAtom());
        }
    }

    private static final long serialVersionUID = 2L;
}
