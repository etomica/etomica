package etomica.species;

import etomica.atom.*;
import etomica.config.IConformation;
import etomica.molecule.IMolecule;
import etomica.molecule.Molecule;
import etomica.space.Space;

import java.util.*;
import java.util.stream.IntStream;

public class SpeciesGeneral implements ISpecies {
    private final AtomType[] atomTypes;
    private final AtomType[] uniqueAtomTypes;
    private final IConformation conformation;
    private final String[] atomNames;
    private final int[][] orientation;
    private final boolean isDynamic;
    private final Space space;
    private int index;

    public SpeciesGeneral(AtomType[] atomTypes,
                          IConformation conformation,
                          String[] atomNames,
                          int[][] orientation,
                          boolean isDynamic,
                          Space space
                          ) {
        this.atomTypes = atomTypes;
        this.conformation = conformation;
        this.atomNames = atomNames;
        this.orientation = orientation;
        this.index = -1;
        LinkedHashSet<AtomType> unique = new LinkedHashSet<>(Arrays.asList(atomTypes));
        this.uniqueAtomTypes = unique.toArray(new AtomType[0]);
        this.isDynamic = isDynamic;
        this.space = space;
    }

    @Override
    public int getIndex() {
        return this.index;
    }

    @Override
    public void setIndex(int newIndex) {
        this.index = newIndex;
    }

    @Override
    public IMolecule makeMolecule() {
        Molecule molecule = new Molecule(this, this.atomTypes.length);
        for (AtomType atomType : this.atomTypes) {
            Atom atom;
            if (atomType instanceof AtomTypeOriented) {
                atom = this.isDynamic ? new AtomOrientedDynamic(this.space, atomType)
                        : new AtomOriented(this.space, atomType);
            } else {
                atom = this.isDynamic ? new AtomLeafDynamic(this.space, atomType)
                        : new Atom(this.space, atomType);
            }
            molecule.addChildAtom(atom);
        }
        this.conformation.initializePositions(molecule.getChildList());
        return molecule;
    }

    @Override
    public int getAtomTypeCount() {
        return this.uniqueAtomTypes.length;
    }

    @Override
    public AtomType getAtomType(int index) {
        return this.uniqueAtomTypes[index];
    }

    @Override
    public List<AtomType> getAtomTypes() {
        return Arrays.asList(this.uniqueAtomTypes);
    }

    @Override
    public void initializeConformation(IMolecule molecule) {
        this.conformation.initializePositions(molecule.getChildList());
    }

    public int getByName(String atomName) {
        for (int i = 0; i < this.atomNames.length; i++) {
            if (this.atomNames[i].equals(atomName) && !atomName.isEmpty()) {
                return i;
            }
        }
        throw new IllegalArgumentException("Name not found");
    }

    /**
     * Gets the index within the molecule of the nth atom of given element string.
     *
     * @param name String matching the Element symbol for the AtomType
     * @param number The nth atom with that type to get, starting at 1.
     * @return the index of that atom within a molecule of this species.
     */
    public int getAtomByTypeName(String name, int number) {
        return IntStream.range(0, this.atomTypes.length)
                .filter(i -> this.atomTypes[i].getElement().getSymbol().equals(name))
                .skip(number - 1)
                .findFirst().orElseThrow(NoSuchElementException::new);
    }

    public int getAtomByTypeName(String name) {
        return this.getAtomByTypeName(name, 1);
    }

    public AtomType getTypeByName(String typeName) {
        for (AtomType type : this.uniqueAtomTypes) {
            if (type.getElement().getSymbol().equals(typeName)) {
                return type;
            }
        }
        throw new IllegalArgumentException("Type not found: " + typeName);
    }

    @Override
    public String toString() {
        return new StringJoiner(", ", SpeciesGeneral.class.getSimpleName() + "[", "]")
                .add("atoms=" + Arrays.toString(atomTypes))
//                .add("uniqueAtomTypes=" + Arrays.toString(uniqueAtomTypes))
//                .add("conformation=" + conformation)
//                .add("atomNames=" + Arrays.toString(atomNames))
//                .add("orientation=" + Arrays.toString(orientation))
//                .add("isDynamic=" + isDynamic)
//                .add("space=" + space)
//                .add("index=" + index)
                .toString();
    }
}
