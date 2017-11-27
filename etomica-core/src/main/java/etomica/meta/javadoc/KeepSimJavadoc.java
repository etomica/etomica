package etomica.meta.javadoc;

import com.github.therapi.runtimejavadoc.RetainJavadoc;

import java.lang.annotation.ElementType;
import java.lang.annotation.Inherited;
import java.lang.annotation.Target;

@RetainJavadoc
@Inherited
@Target(ElementType.TYPE)
public @interface KeepSimJavadoc {

}


