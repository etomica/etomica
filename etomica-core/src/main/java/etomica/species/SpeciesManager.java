package etomica.species;

import etomica.atom.AtomType;
import etomica.chem.elements.IElement;

import java.util.*;

public final class SpeciesManager {

    private final ISpecies[] speciesArray;
    private final List<ISpecies> speciesList;
    private final int atomTypeCount;
    private final boolean isPureAtoms;

    // TODO should it be public? would need to handle setIndex
    private SpeciesManager(ISpecies[] speciesArray) {
        this.speciesArray = speciesArray;
        this.speciesList = Arrays.asList(speciesArray);
        this.atomTypeCount = Arrays.stream(speciesArray).mapToInt(ISpecies::getUniqueAtomTypeCount).sum();
        this.isPureAtoms = Arrays.stream(speciesArray).allMatch(s -> s.getLeafAtomCount() == 1);
    }

    public static Builder builder() {
        return new Builder();
    }

    public ISpecies[] getSpeciesArray() {
        return this.speciesArray;
    }

    public List<ISpecies> getSpeciesList() {
        return this.speciesList;
    }

    public ISpecies getSpecies(int idx) {
        return speciesArray[idx];
    }

    public int getSpeciesCount() {
        return this.speciesArray.length;
    }

    public int getAtomTypeCount() {
        return this.atomTypeCount;
    }

    public boolean isPureAtoms() {
        return this.isPureAtoms;
    }

    public static class Builder {
        private final Map<String, IElement> elementSymbols = new HashMap<>();
        private final List<ISpecies> speciesList = new ArrayList<>();
        private int atomTypeCount = 0;

        public SpeciesManager build() {
            return new SpeciesManager(speciesList.toArray(new ISpecies[0]));
        }

        public Builder addSpecies(ISpecies species) {
            Objects.requireNonNull(species);
            if (speciesList.contains(species)) {
                throw new IllegalArgumentException("Species is already added: " + species);
            }

            species.setIndex(speciesList.size());
            speciesList.add(species);

            for (AtomType atomType : species.getUniqueAtomTypes()) {
                atomType.setIndex(atomTypeCount);
                atomTypeCount++;

                IElement newElement = atomType.getElement();
                elementSymbols.merge(newElement.getSymbol(), newElement, (oldElement, newEl) -> {
                    if (oldElement != newEl) {
                        throw new IllegalStateException(
                                "Element symbol " + newEl.getSymbol() + " already exists in this simulation as a different element"
                        );
                    }
                    return newEl;
                });
            }

            return this;
        }

        public Builder addSpecies(ISpecies... species) {
            Arrays.stream(species).forEach(this::addSpecies);
            return this;
        }

        /**
         * Method to allow generation of unique string to identify Elements in the Simulation.
         *
         * @param symbolBase the base string.
         * @return an Element symbol starting with symbolBase that does not yet
         * exist in the Simulation.  Return values will be like "base0, base1, base2..."
         */
        public String makeUniqueElementSymbol(String symbolBase) {
            int n = 0;
            while (elementSymbols.containsKey(symbolBase + n)) {
                n++;
            }
            // reserve this symbol so future calls to makeUniqueElementSymbol won't return it
            // this will get replaced by the actual Element when it gets added
            elementSymbols.put(symbolBase + n, null);
            return symbolBase + n;
        }

    }
}
