package etomica.GasMOP;

import etomica.atom.AtomType;
import etomica.chem.elements.*;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

public class AtomTypeFactory {
    // cache: "C_1" -> AtomType instance
    // cache: "C_1" -> AtomType
    private static final Map<String, AtomType> CACHE = new ConcurrentHashMap<>();

    // element symbol -> element singleton
    private static final Map<String, IElement> ELEMENTS;

    static {
        Map<String, IElement> m = new HashMap<>();
        m.put("C", Carbon.INSTANCE);
        m.put("O", Oxygen.INSTANCE);
        m.put("H", Hydrogen.INSTANCE);
        m.put("N", Nitrogen.INSTANCE);
        // add more if needed
        ELEMENTS = Collections.unmodifiableMap(m);
    }

    public AtomTypeFactory() {}

    /** Returns a cached AtomType for labels like "C_5", "O_2" */
    public AtomType fromLabel(String label) {
        AtomType at = CACHE.get(label);
        if (at != null) return at;

        // create lazily
        at = createAtomType(label);
        CACHE.put(label, at);
        return at;
    }

    public AtomType createAtomType(String label) {
        int idx = label.indexOf('_');
        String symbol = (idx >= 0) ? label.substring(0, idx) : label;

        IElement elem = ELEMENTS.get(symbol);
        if (elem == null) {
            throw new RuntimeException("Unknown element symbol in atom type: " + label);
        }

        return new AtomType(elem, label);
    }

    public static void main(String[] args) {
    }
}
