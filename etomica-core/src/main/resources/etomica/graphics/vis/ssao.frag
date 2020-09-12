#version 330 core
out float FragColor;

in vec2 TexCoords;

uniform sampler2D gPosition;
uniform sampler2D gNormal;
uniform sampler2D texNoise;

uniform vec3 ssaoKernel[64];
uniform mat4 projection;
uniform float width;
uniform float height;

const float radius = 0.5;
const float bias = 0.025;

void main() {
    FragColor = 1;
    vec3 normal = texture(gNormal, TexCoords).rgb;
    if (normal == vec3(0)) { return; } // Don't know if this is a good idea or not
    vec2 noiseScale = vec2(width / 4.0, height / 4.0);
    vec3 fragPos = texture(gPosition, TexCoords).xyz;
    vec3 randomVec = texture(texNoise, TexCoords * noiseScale).rgb;

    vec3 tangent = normalize(randomVec - normal * dot(randomVec, normal));
    vec3 bitangent = cross(normal, tangent);
    mat3 TBN = mat3(tangent, bitangent, normal);

    float occlusion = 0.0;
    for (int i = 0; i < 32; i++) {
        vec3 occSample = TBN * ssaoKernel[i];
        occSample = fragPos + occSample * radius;

        vec4 offset = vec4(occSample, 1.0);
        offset = projection * offset; // view to clip-space
        offset.xyz /= offset.w; // perspective divide
        offset.xyz = offset.xyz * 0.5 + 0.5; // transform to range 0.0-1.0
        float sampleDepth = texture(gPosition, offset.xy).z;
        float rangeCheck = smoothstep(0.0, 1.0, radius / abs(fragPos.z - sampleDepth));
        occlusion += (sampleDepth >= occSample.z + bias ? 1.0 : 0.0) * rangeCheck;
    }
    occlusion = 1.0 - (occlusion / 32);
    const float power = 2.0;
    FragColor = pow(occlusion, power);
//    FragColor=1;
}
