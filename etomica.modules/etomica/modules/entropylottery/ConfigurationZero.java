package etomica.modules.entropylottery;

import etomica.action.AtomActionTranslateTo;
import etomica.atom.IAtom;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.config.Configuration;
import etomica.phase.Phase;
import etomica.space.IVector;

/**
 * Configuration that simply puts all the Atoms at 0, or -0.5 if the box size
 * is even (which makes the entropy lottery binning happy).
 * @author Andrew Schultz
 */
public class ConfigurationZero extends Configuration {


    public ConfigurationZero() {
        super();
    }

    public void initializeCoordinates(Phase phase) {
        AtomActionTranslateTo atomActionTranslateTo = new AtomActionTranslateTo(phase.getSpace());
        IVector work = phase.getSpace().makeVector();
        work.E(0.0);
        int intD = (int)Math.round(phase.getBoundary().getDimensions().x(0));
        if (intD % 2 == 0) {
            work.E(-0.5);
        }
        atomActionTranslateTo.setDestination(work);

        AtomIteratorAllMolecules iterator = new AtomIteratorAllMolecules(phase);
        iterator.reset();
        for (IAtom a = iterator.nextAtom(); a != null; a = iterator.nextAtom()) {
           atomActionTranslateTo.actionPerformed(iterator.nextAtom());
        }
    }

    private static final long serialVersionUID = 2L;
}
