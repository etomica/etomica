package etomica.server.representations;

import com.fasterxml.jackson.annotation.JsonProperty;

public class SimClassInfo {
    private final @JsonProperty String className;
    private final @JsonProperty String javadoc;

    public SimClassInfo(String className, String javadoc) {
        this.className = className;
        this.javadoc = javadoc;
    }
}
