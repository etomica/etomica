package etomica.species;

import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.IConformation;
import etomica.molecule.IMolecule;
import etomica.space.Space;
import etomica.space.Vector;

import java.lang.reflect.Field;
import java.util.List;

public class SpeciesGeneralMutable implements ISpecies {

    private SpeciesGeneral species;
    private int index = -1;

    public SpeciesGeneralMutable(SpeciesGeneral species) {
        this.species = species;
        setParent(this.species);
    }

    public void setSpecies(SpeciesGeneral species) {
        this.species = species;
        setParent(this.species);
    }

    private void setParent(SpeciesGeneral species) {
        species.getUniqueAtomTypes().forEach(atomType -> atomType.setSpecies(this));
    }

    public SpeciesGeneral getInnerSpecies() {
        return this.species;
    }

    @Override
    public int getIndex() {
        return index;
    }

    @Override
    public void setIndex(int newIndex) {
        this.index = newIndex;
    }

    @Override
    public IMolecule makeMolecule() {
        IMolecule mol = species.makeMolecule();
        try {
            Field species = mol.getClass().getDeclaredField("species");
            species.setAccessible(true);
            species.set(mol, this);
        } catch (NoSuchFieldException | IllegalAccessException e) {
            throw new RuntimeException(e);
        }
        return mol;
    }

    @Override
    public IMolecule initMolecule(Box box, int molIdx, int atomIdxStart) {
        IMolecule mol = species.initMolecule(box, molIdx, atomIdxStart);
        try {
            Field species = mol.getClass().getDeclaredField("species");
            species.setAccessible(true);
            species.set(mol, this);
        } catch (NoSuchFieldException | IllegalAccessException e) {
            throw new RuntimeException(e);
        }
        return mol;
    }

    @Override
    public int getUniqueAtomTypeCount() {
        return species.getUniqueAtomTypeCount();
    }

    @Override
    public int getLeafAtomCount() {
        return species.getLeafAtomCount();
    }

    @Override
    public AtomType getAtomType(int index) {
        return species.getAtomType(index);
    }

    @Override
    public AtomType getLeafType() {
        return species.getLeafType();
    }

    @Override
    public List<AtomType> getUniqueAtomTypes() {
        return species.getUniqueAtomTypes();
    }

    @Override
    public List<AtomType> getAtomTypes() {
        return species.getAtomTypes();
    }

    @Override
    public void initializeConformation(IMolecule molecule) {
        species.initializeConformation(molecule);
    }

    @Override
    public IConformation getConformation() {
        return species.getConformation();
    }

    @Override
    public int getByName(String atomName) {
        return species.getByName(atomName);
    }

    @Override
    public int getAtomByTypeName(String name, int number) {
        return species.getAtomByTypeName(name, number);
    }

    @Override
    public int getAtomByTypeName(String name) {
        return species.getAtomByTypeName(name);
    }

    @Override
    public AtomType getTypeByName(String typeName) {
        return species.getTypeByName(typeName);
    }

    @Override
    public Vector getMomentOfInertia() {
        return species.getMomentOfInertia();
    }

    @Override
    public double getMass() {
        return species.getMass();
    }

    @Override
    public String toString() {
        return species.toString();
    }

    @Override
    public boolean isDynamic() {
        return species.isDynamic();
    }

    @Override
    public Space getSpace() {
        return species.getSpace();
    }
}
