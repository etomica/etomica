#version 330 core

out vec4 FragColor;

in vec2 TexCoords;

uniform sampler2D gPosition;
uniform sampler2D gNormal;
uniform sampler2D gAlbedoSpec;
uniform sampler2D ssao;

struct DirLight {
    vec3 direction;

    vec3 ambient;
    vec3 diffuse;
    vec3 specular;
};
uniform DirLight light;
uniform vec3 backgroundColor;

void main() {

    vec3 frag_pos = texture(gPosition, TexCoords).rgb;
    vec3 normal = texture(gNormal, TexCoords).rgb;
    vec3 diffuse = texture(gAlbedoSpec, TexCoords).rgb;
    float specular = texture(gAlbedoSpec, TexCoords).a;
    float ao = texture(ssao, TexCoords).r;

    vec3 viewDir = normalize(frag_pos);

    vec3 lightDir = normalize(-light.direction);
    vec3 halfwayDir = normalize(lightDir + viewDir);

    float diffuseFac = max(dot(normal, lightDir), 0.0);

    vec3 reflectDir = reflect(-lightDir, normal);
    float specularFac = pow(max(dot(normal, halfwayDir), 0.0), 16.0);

    vec3 ambientColor = light.ambient * diffuse * ao;
    vec3 diffuseColor = light.diffuse * diffuseFac * diffuse;
    vec3 spec = light.specular * specularFac * specular;

    float edgeFac = dot(-viewDir, normal);
    float fac = smoothstep(0.3, 0.4, edgeFac);
    vec3 color = (ambientColor + diffuseColor + spec);
    FragColor = vec4(color * fac, 1.0);
    if (normal == vec3(0)) {
        FragColor = vec4(backgroundColor, 0);
    }
}
