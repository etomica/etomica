package etomica.units.dimensions;

import etomica.units.dimensions.Dimension;

import java.lang.annotation.ElementType;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;

@Retention(RetentionPolicy.RUNTIME)
@Target(ElementType.METHOD)
public @interface Dimensioned {
    Class<? extends Dimension> dimension();
}
