package etomica.species;

import etomica.atom.*;
import etomica.box.Box;
import etomica.box.storage.Tokens;
import etomica.config.IConformation;
import etomica.molecule.*;
import etomica.space.IOrientation;
import etomica.space.Space;
import etomica.space.Vector;

import java.util.*;
import java.util.stream.IntStream;

public class SpeciesGeneral implements ISpecies {
    private final AtomType[] atomTypes;
    private final AtomType[] uniqueAtomTypes;
    private final IConformation conformation;
    private final String[] atomNames;
    private final int[][] orientation;
    private final boolean isDynamic;
    private final boolean isMoleculeOriented;
    private final Space space;
    private final Vector momentOfInertia;
    private final double mass;
    private final AtomFactory atomFactory;
    private int index;

    public SpeciesGeneral(AtomType[] atomTypes,
                          IConformation conformation,
                          String[] atomNames,
                          int[][] orientation,
                          boolean isDynamic,
                          boolean isMoleculeOriented,
                          Vector momentOfInertia,
                          Space space,
                          AtomFactory atomFactory) {
        this.atomTypes = atomTypes;
        this.conformation = conformation;
        this.atomNames = atomNames;
        this.orientation = orientation;
        this.isMoleculeOriented = isMoleculeOriented;
        this.atomFactory = atomFactory;
        this.index = -1;
        LinkedHashSet<AtomType> unique = new LinkedHashSet<>(Arrays.asList(atomTypes));
        this.uniqueAtomTypes = unique.toArray(new AtomType[0]);
        this.isDynamic = isDynamic;
        this.space = space;

        for (AtomType atomType : unique) {
            atomType.setSpecies(this);
        }

        this.mass = Arrays.stream(this.atomTypes).mapToDouble(atomType -> atomType.getMass()).sum();
        if (momentOfInertia == null) {
            this.momentOfInertia = calcMomentOfInertia();
        } else {
            this.momentOfInertia = momentOfInertia;
        }
    }

    public static SpeciesGeneral monatomic(Space space, AtomType atomType) {
        return SpeciesGeneral.monatomic(space, atomType, false);
    }

    public static SpeciesGeneral monatomic(Space space, AtomType atomType, boolean isDynamic) {
        return new SpeciesBuilder(space)
                .addAtom(atomType, space.makeVector())
                .setDynamic(isDynamic)
                .build();
    }

    private Vector calcMomentOfInertia() {
        if (space.D() != 3) {
            return null;
        }
        // make a pretend molecule and calculate its moment of inertia
        IAtomList dummyAtoms = new AtomArrayList(this.getLeafAtomCount());
        for (AtomType atomType : this.atomTypes) {
            dummyAtoms.add(new Atom(space, atomType, space.makeVector()));
        }
        conformation.initializePositions(dummyAtoms);
        Vector com = MoleculePositionCOM.centerOfMass(dummyAtoms, space.makeVector());
        double[] I = new double[3];
        Vector xWork = space.makeVector();
        Vector moment = space.makeVector();
        for (IAtom atom : dummyAtoms) {
            xWork.Ev1Mv2(atom.getPosition(), com);
            double atomMass = atom.getType().getMass();
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    if (j == k) continue;
                    I[j] += atomMass * xWork.getX(k) * xWork.getX(k);
                }
            }
        }
        moment.E(I);
        return moment;
    }

    @Override
    public int getIndex() {
        return this.index;
    }

    @Override
    public void setIndex(int newIndex) {
        this.index = newIndex;
    }

    public IMolecule initMolecule(Box box, int molIdx, int atomIdxStart) {
        Molecule molecule;
        if (this.isMoleculeOriented) {
            IOrientation orientation = box.getMolOrientations(Tokens.ORIENTATION_FULL).create(molIdx);
            Vector position = box.getMolVectors(Tokens.POSITION).create(molIdx);
            if (this.isDynamic) {
                Vector angularMomentum = box.getMolVectors(Tokens.ANGULAR_MOMENTUM).create(molIdx);
                Vector velocity = box.getMolVectors(Tokens.VELOCITY).create(molIdx);
                molecule = new MoleculeOrientedDynamic(this, orientation, position, angularMomentum, velocity);
            } else {
                molecule = new MoleculeOriented(this, orientation, position);
            }
        } else {
            molecule = new Molecule(this);
        }
        molecule.setGlobalIndex(molIdx);

        for (int i = 0; i < atomTypes.length; i++) {
            AtomType type = atomTypes[i];
            int atomIdx = atomIdxStart + i;
            IAtom atom = this.atomFactory.makeAtom(type, box, atomIdx);
            molecule.addChildAtom(atom);
        }
        this.conformation.initializePositions(molecule.getChildList());
        return molecule;
    }

    @Override
    public boolean isDynamic() {
        return isDynamic;
    }

    @Override
    public Space getSpace() {
        return space;
    }

    @Override
    public int getUniqueAtomTypeCount() {
        return this.uniqueAtomTypes.length;
    }

    @Override
    public int getLeafAtomCount() {
        return this.atomTypes.length;
    }

    @Override
    public AtomType getAtomType(int index) {
        return this.uniqueAtomTypes[index];
    }

    @Override
    public AtomType getLeafType() {
        if (this.uniqueAtomTypes.length > 1) {
            throw new RuntimeException("Species " + this.toString() + " has more than one atom type");
        }
        return getAtomType(0);
    }

    @Override
    public List<AtomType> getUniqueAtomTypes() {
        return Arrays.asList(this.uniqueAtomTypes);
    }

    @Override
    public List<AtomType> getAtomTypes() {
        return Arrays.asList(this.atomTypes);
    }

    @Override
    public void initializeConformation(IMolecule molecule) {
        this.conformation.initializePositions(molecule.getChildList());
    }

    public IConformation getConformation() {
        return this.conformation;
    }

    @Override
    public int getByName(String atomName) {
        for (int i = 0; i < this.atomNames.length; i++) {
            if (this.atomNames[i].equals(atomName) && !atomName.isEmpty()) {
                return i;
            }
        }
        throw new IllegalArgumentException("Name not found");
    }

    @Override
    public int getAtomByTypeName(String name, int number) {
        return IntStream.range(0, this.atomTypes.length)
                .filter(i -> this.atomTypes[i].getName().equals(name))
                .skip(number - 1)
                .findFirst().orElseThrow(NoSuchElementException::new);
    }

    @Override
    public int getAtomByTypeName(String name) {
        return this.getAtomByTypeName(name, 1);
    }

    @Override
    public AtomType getTypeByName(String typeName) {
        for (AtomType type : this.uniqueAtomTypes) {
            if (type.getName().equals(typeName)) {
                return type;
            }
        }
        throw new IllegalArgumentException("Type not found: " + typeName);
    }

    @Override
    public Vector getMomentOfInertia() {
        return this.momentOfInertia;
    }

    @Override
    public double getMass() {
        return mass;
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

    public static AtomFactory defaultAtomFactory(Space space, boolean isDynamic) {
        return (atomType, box, id) -> {
            Vector position = box.getAtomVectors(Tokens.POSITION).create(id);
            if (atomType instanceof AtomTypeOriented) {
                IOrientation orientation = box.getAtomOrientations(space.D() == 3 ? Tokens.ORIENTATION_FULL : Tokens.ORIENTATION).create(id);
                if (isDynamic) {
                    Vector velocity = box.getAtomVectors(Tokens.VELOCITY).create(id);
                    Vector angularVel = box.getAtomVectors(Tokens.ANGULAR_VELOCITY).create(id);
                    return new AtomOrientedDynamic(space, atomType, orientation, position, velocity, angularVel);
                } else {
                    return new AtomOriented(space, atomType, position, orientation);
                }
            } else {
                if (isDynamic) {
                    Vector velocity = box.getAtomVectors(Tokens.VELOCITY).create(id);
                    return new AtomLeafDynamic(space, atomType, position, velocity);
                } else {
                    return new Atom(space, atomType, position);
                }
            }
        };
    }

    public interface AtomFactory {
        IAtom makeAtom(AtomType atomType, Box box, int id);
    }
}
