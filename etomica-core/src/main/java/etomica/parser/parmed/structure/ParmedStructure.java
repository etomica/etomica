package etomica.parser.parmed.structure;

import com.fasterxml.jackson.annotation.*;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.PropertyNamingStrategy;
import com.fasterxml.jackson.databind.annotation.JsonNaming;
import etomica.parser.parmed.Views;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

@JsonIgnoreProperties(ignoreUnknown = false)
@JsonNaming(PropertyNamingStrategy.SnakeCaseStrategy.class)
public class ParmedStructure {

    public String title;

    @JsonProperty("_combining_rule")
    public String combiningRule;

    @JsonProperty("_coordinates")
    public double[][][] coordinates;

    @JsonProperty("_box")
    public double[][] box;

    public List<Residue> residues;

    public List<ParameterType> bondTypes;
    public List<ParameterType> angleTypes;
    public List<ParameterType> adjustTypes;
    public List<ParameterType> rbTorsionTypes;

    public int[][] bonds;
    public int[][] angles;
    public int[][] adjusts;
    public int[][] rbTorsions;


    public Map<String, Object> defaults;

    @JsonAnySetter
    public final Map<String, Object> otherProperties = new HashMap<>();

    @JsonNaming(PropertyNamingStrategy.SnakeCaseStrategy.class)
    public static class Atom {
        public String name;
        public String type;
        public double mass;
        public double xx;
        public double xy;
        public double xz;
        public int atomicNumber;

        @JsonView(Views.AtomTypes.class)
        public AtomType atomType;

        @JsonProperty("_charge")
        public double charge;


        @JsonAnySetter
        public final Map<String, Object> properties = new HashMap<>();
    }

    @JsonNaming(PropertyNamingStrategy.SnakeCaseStrategy.class)
    public static class AtomType {
        public String name;
        public double mass;
        public Double number;
        public int atomicNumber;


        @JsonAnySetter
        public final Map<String, Object> properties = new HashMap<>();
    }

    @JsonNaming(PropertyNamingStrategy.SnakeCaseStrategy.class)
    public static class Residue {
        public String name;
        public int number;
        public String chain;
        public List<Atom> atoms;

        @JsonAnySetter
        public final Map<String, Object> properties = new HashMap<>();
    }

    public static void main(String[] args) throws IOException {
        ObjectMapper om = new ObjectMapper();
        ParmedStructure p = om
                .readerWithView(Views.Default.class)
                .forType(ParmedStructure.class)
                .readValue(new File("/home/alex/workspace/parmed_json/pdb.json"));
        System.out.println();
    }
}


