package etomica.server.representations;

import com.fasterxml.jackson.annotation.JsonProperty;

import java.util.List;
import java.util.Map;

public class MeterInfo {
    private final @JsonProperty String className;
    private final @JsonProperty List<String> constructorParamTypes;
    private final @JsonProperty Map<String, List<Long>> classOptions;
    private final @JsonProperty Map<String, String> propertyTypes;

    public MeterInfo(String className, List<String> constructorParamTypes, Map<String, List<Long>> classOptions, Map<String, String> propertyTypes) {
        this.className = className;
        this.constructorParamTypes = constructorParamTypes;
        this.classOptions = classOptions;
        this.propertyTypes = propertyTypes;
    }
}
