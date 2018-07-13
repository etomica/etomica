module etomica.core {
    requires jackson.annotations;
    requires java.desktop;
    requires jackson.databind;
    requires therapi.runtime.javadoc;
    requires etomica.graphics3D;
    requires jama;
    exports etomica.action;
    exports etomica.atom;
    exports etomica.data;
    exports etomica.data.history;
    exports etomica.data.meter;
    exports etomica.data.types;
    exports etomica.graphics;
    exports etomica.math.geometry;
    exports etomica.modifier;
    exports etomica.units;
    exports etomica.space3d;
    exports etomica.space;
    exports etomica.units.dimensions;
}
