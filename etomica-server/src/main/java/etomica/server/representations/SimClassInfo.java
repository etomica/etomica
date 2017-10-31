package etomica.server.representations;

import com.fasterxml.jackson.annotation.JsonProperty;
import com.github.therapi.runtimejavadoc.ClassJavadoc;
import com.github.therapi.runtimejavadoc.Comment;
import com.github.therapi.runtimejavadoc.RuntimeJavadoc;
import org.jsoup.Jsoup;
import org.jsoup.safety.Whitelist;

import java.util.Optional;

public class SimClassInfo {
    private final @JsonProperty String className;
    private final @JsonProperty String javadoc;

    public SimClassInfo(String className, String javadoc) {
        this.className = className;
        this.javadoc = javadoc;
    }

    public static SimClassInfo forClass(Class<?> cls) {

        Optional<ClassJavadoc> javadoc = RuntimeJavadoc.getJavadoc(cls);
        String comment = javadoc.map(ClassJavadoc::getComment).map(Comment::toString).orElse("");
        String sanitizedComment = Jsoup.clean(comment, Whitelist.basic());
        return new SimClassInfo(cls.getCanonicalName(), sanitizedComment);
    }
}
