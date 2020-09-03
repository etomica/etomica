package etomica.species;

import etomica.atom.AtomType;
import etomica.molecule.IMolecule;
import etomica.species.ISpecies;
import etomica.species.SpeciesGeneral;

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
        species.getAtomTypes().forEach(atomType -> atomType.setSpecies(this));
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
    public int getAtomTypeCount() {
        return species.getAtomTypeCount();
    }

    @Override
    public AtomType getAtomType(int index) {
        return species.getAtomType(index);
    }

    @Override
    public List<AtomType> getAtomTypes() {
        return species.getAtomTypes();
    }

    @Override
    public void initializeConformation(IMolecule molecule) {
        species.initializeConformation(molecule);
    }
}
