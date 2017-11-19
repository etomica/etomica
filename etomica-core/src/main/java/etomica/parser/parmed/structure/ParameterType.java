package etomica.parser.parmed.structure;

import com.fasterxml.jackson.annotation.JsonAnySetter;
import com.fasterxml.jackson.annotation.JsonProperty;

import java.util.HashMap;
import java.util.Map;

public class ParameterType {
    @JsonProperty
    public boolean used;

    @JsonProperty
    public Object penalty;

    @JsonProperty
    public int _idx;

    public final Map<String, Double> properties = new HashMap<>();

    @JsonAnySetter
    public void setProperty(String name, Double prop) {
        this.properties.put(name, prop);
    }
}
