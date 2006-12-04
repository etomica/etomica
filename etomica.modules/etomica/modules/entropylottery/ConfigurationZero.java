package etomica.modules.entropylottery;

import etomica.action.AtomActionTranslateTo;
import etomica.atom.AtomArrayList;
import etomica.atom.iterator.AtomIteratorArrayListCompound;
import etomica.config.Configuration;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * Configuration that simply puts all the Atoms at 0, or -0.5 if the box size
 * is even (which makes the entropy lottery binning happy).
 * @author Andrew Schultz
 */
public class ConfigurationZero extends Configuration {


    public ConfigurationZero(Space space) {
        super(space);
    }

    protected void initializePositions(AtomArrayList[] atomList) {
        AtomActionTranslateTo atomActionTranslateTo = new AtomActionTranslateTo(space);
        Vector work = space.makeVector();
        work.E(0.0);
        int intD = (int)Math.round(dimensions[0]);
        if (intD % 2 == 0) {
            work.E(-0.5);
        }
        atomActionTranslateTo.setDestination(work);

        AtomIteratorArrayListCompound iterator = new AtomIteratorArrayListCompound(atomList);
        iterator.reset();
        while (iterator.hasNext()) {
           atomActionTranslateTo.actionPerformed(iterator.nextAtom());
        }
    }

    private static final long serialVersionUID = 1L;
}
