package etomica.species;

import etomica.atom.AtomType;
import etomica.config.ConformationGeneric;
import etomica.config.IConformation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

public class SpeciesBuilder {
    private final List<AtomType> atomTypes;
    private final List<String> atomNames;
    private final List<Vector> atomPositions;
    private IConformation conformation;
    private final List<int[]> orientation;
    private boolean isDynamic;
    private final Space space;

    public SpeciesBuilder(Space space) {
        this.atomTypes = new ArrayList<>();
        this.atomPositions = new ArrayList<>();
        this.conformation = null;
        this.orientation = new ArrayList<>(2);
        this.atomNames = new ArrayList<>();
        this.isDynamic = false;
        this.space = space;
    }

    public SpeciesBuilder setDynamic(boolean isDynamic) {
        this.isDynamic = isDynamic;
        return this;
    }

    public SpeciesBuilder withConformation(IConformation conformation) {
        this.conformation = conformation;
        this.check();
        return this;
    }

    public SpeciesBuilder addAtom(AtomType type, Vector position) {
        return this.addAtom(type, position, "");
    }

    public SpeciesBuilder addAtom(AtomType type, Vector position, String name) {
        if (!name.isEmpty() && new HashSet<>(this.atomNames).contains(name)) {
            throw new IllegalArgumentException("Atom name already used");
        }
        this.atomTypes.add(type);
        this.atomPositions.add(position);
        this.atomNames.add(name);
        this.check();
        return this;
    }

    public SpeciesBuilder addCount(AtomType type, int count) {
        for (int i = 0; i < count; i++) {
            this.atomTypes.add(type);
            this.atomNames.add("");
        }
        this.check();
        return this;
    }

    public SpeciesBuilder withDirection(int fromAtom, int toAtom) {
        if (fromAtom == toAtom) {
            throw new IllegalStateException("Direction vector cannot be zero");
        }
        if (orientation.size() == 2) {
            throw new IllegalStateException("Only two directions allowed");
        }
        this.orientation.add(new int[] {fromAtom, toAtom});
        return this;
    }

    public SpeciesGeneral build() {
        // TODO check conditions
        if (this.conformation == null) {
            this.conformation = new ConformationGeneric(this.atomPositions.toArray(new Vector[0]));
        }

        return new SpeciesGeneral(
                this.atomTypes.toArray(new AtomType[0]),
                this.conformation,
                this.atomNames.toArray(new String[0]),
                this.orientation.toArray(new int[0][]),
                this.isDynamic,
                this.space
        );
    }

    private void check() {
        if (!this.atomPositions.isEmpty() && this.conformation != null) {
            throw new IllegalStateException("");
        }

        if (this.atomPositions.size() != this.atomTypes.size() && !atomPositions.isEmpty()) {
            throw new IllegalStateException("Cannot specify an atom count and add single atoms");
        }
    }

}

