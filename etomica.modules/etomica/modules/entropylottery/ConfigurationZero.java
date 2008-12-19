package etomica.modules.entropylottery;

import etomica.action.MoleculeActionTranslateTo;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IVectorMutable;
import etomica.atom.iterator.MoleculeIteratorAllMolecules;
import etomica.config.Configuration;
import etomica.space.ISpace;

/**
 * Configuration that simply puts all the Atoms at 0, or -0.5 if the box size
 * is even (which makes the entropy lottery binning happy).
 * @author Andrew Schultz
 */
public class ConfigurationZero implements Configuration, java.io.Serializable {

	private final ISpace space;

    public ConfigurationZero(ISpace _space) {
        super();
        this.space = _space;
    }

    public void initializeCoordinates(IBox box) {
        MoleculeActionTranslateTo atomActionTranslateTo = new MoleculeActionTranslateTo(space);
        IVectorMutable work = space.makeVector();
        work.E(0.0);
        int intD = (int)Math.round(box.getBoundary().getDimensions().x(0));
        if (intD % 2 == 0) {
            work.E(-0.5);
        }
        atomActionTranslateTo.setDestination(work);

        MoleculeIteratorAllMolecules iterator = new MoleculeIteratorAllMolecules(box);
        iterator.reset();
        for (IMolecule a = iterator.nextMolecule(); a != null; a = iterator.nextMolecule()) {
           atomActionTranslateTo.actionPerformed(iterator.nextMolecule());
        }
    }

    private static final long serialVersionUID = 2L;
}
