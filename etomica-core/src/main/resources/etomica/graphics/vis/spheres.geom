#version 330 core

/*
    From "2D Polyhedral Bounds of a Clipped, Perspective-Projected 3D Sphere"
    http://jcgt.org/published/0002/02/05/paper.pdf
*/

/** The number of sides in the bounding polygon. Must be even. */
#define N 4
/** The number of axes to caclulate bounds along. Two bounds per axis */
#define N_AXES (N/2)
#define PI (3.1415926)

layout (points) in;
layout (triangle_strip, max_vertices = N) out;

in float vRadius[];

uniform mat4 projection;
uniform float nearZ;

flat out vec3 sphere_center;
flat out float radius_sq;
smooth out vec3 frag_pos;


/** 2D-line from point and direction */
struct line2D {
    vec2 point;
    vec2 direction;
};

/** Simple line-line intersection in 2D. We don't handle
    parallel lines, since they cannot arise in our implementation */
vec2 intersect(line2D L1, line2D L2){
    float denominator = (L1.direction.x * L2.direction.y) -
    (L1.direction.y * L2.direction.x);
    float leftTerm  =   (L1.point.x + L1.direction.x) * L1.point.y -
    (L1.point.y + L1.direction.y) * L1.point.x;
    float rightTerm =   (L2.point.x + L2.direction.x) * L2.point.y -
    (L2.point.y + L2.direction.y) * L2.point.x;
    vec2 numerator = leftTerm * L2.direction - rightTerm * L1.direction;
    return (numerator / denominator);
}


float square(float x){
    return x*x;
}

/**
    Calculates the unprojected upper and lower bounds along axis determined by phi.
    perpendicularDirection is the direction in the xy plane perpendicular to the axis.
    The two lines running in perpendicularDirection through U and L determine two bounding lines
    in a screen-space bounding polygon for the sphere(center,radius) projected onto the screen.
    center, U, and L are all in camera-space.

    Precondition: sphere not culled by near-plane
*/
void getBoundsForPhi(in float phi, in vec3 center, in float radius, in float nearZ, out vec2 perpendicularDirection, out vec3 U, out vec3 L){
    bool trivialAccept = (center.z + radius) < nearZ; // Entirely in back of nearPlane (Trivial Accept)

    vec3 a = vec3(cos(phi), sin(phi), 0);
    perpendicularDirection.x = -a.y;
    perpendicularDirection.y = a.x;

    // given in coordinates (a,z), where a is in the direction of the vector a, and z is in the standard z direction
    vec2 projectedCenter = vec2(dot(a, center), center.z);
    vec2 bounds_az[2];

    float projCenterSqLength = dot(projectedCenter, projectedCenter);
    float tSquared = dot(projectedCenter, projectedCenter) - square(radius);
    float costheta, sintheta;

    if(tSquared >  0) { // Camera is outside sphere
        // Distance to the tangent points of the sphere (points where a vector from the camera are tangent to the sphere) (calculated a-z space)
        float t = sqrt(tSquared);
        float invCLength = inversesqrt(projCenterSqLength);

        // Theta is the angle between the vector from the camera to the center of the sphere and the vectors from the camera to the tangent points
        costheta = t * invCLength;
        sintheta = radius * invCLength;
    }
    float sqrtPart;
    if(!trivialAccept) sqrtPart = sqrt(square(radius) - square(nearZ - projectedCenter.y));

    for( int i = 0; i < 2; ++i ) {
        // Depending on the compiler and platform, it may be possible to optimize for
        // performance by expanding out this caculation and using temporary variables
        // for costheta^2 and costheta*sintheta
        if(tSquared >  0) {
            // Matrices are column-major in GLSL
            mat2 rotateTheta = mat2(    costheta,   sintheta,
            -sintheta,   costheta);
            bounds_az[i] = costheta * (rotateTheta * projectedCenter);
        }

        if(!trivialAccept && (tSquared <= 0 || bounds_az[i].y > nearZ)) {
            bounds_az[i].x = projectedCenter.x + sqrtPart;
            bounds_az[i].y = nearZ;
        }
        sintheta *= -1; // negate theta for B
        sqrtPart *= -1; // negate sqrtPart for B
    }
    U   = bounds_az[0].x * a;
    U.z = bounds_az[0].y;
    L   = bounds_az[1].x * a;
    L.z = bounds_az[1].y;
}

void main() {

    float radius = vRadius[0];
    vec3 C = gl_in[0].gl_Position.xyz;

    if (C.z - radius >= nearZ) { return; }

    // Duplicate the first line into the last spot to avoid modular arithmetic
    line2D boundingLines[N + 1];
    float invNAxes = 1.0 / N_AXES;

    // The plane we draw the bounds on
    float maxZ = min(nearZ, C.z + radius);

    vec3 axesBounds[N];
    for(int i = 0; i < N_AXES; i++) {
        float phi = (i * PI) * invNAxes;
        getBoundsForPhi(phi, C, radius, nearZ, boundingLines[i].direction, axesBounds[i], axesBounds[i + N_AXES]);
        boundingLines[i + N_AXES].direction = boundingLines[i].direction;
    }

    for (int i = 0; i < N; i++) {
        boundingLines[i].point = axesBounds[i].xy * (maxZ / axesBounds[i].z);
    }

    boundingLines[N] = boundingLines[0];

    radius_sq = radius * radius;
    sphere_center = gl_in[0].gl_Position.xyz;

    for (int i = 0; i < N_AXES; i++) {
        int j = N - i - 1;
        vec4 pos = vec4(intersect(boundingLines[j], boundingLines[j + 1]), maxZ, 1.0f);
        frag_pos = pos.xyz;
        gl_Position = projection * pos;
        EmitVertex();
        pos = vec4(intersect(boundingLines[i], boundingLines[i + 1]), maxZ, 1.0f);
        frag_pos = pos.xyz;
        gl_Position = projection * pos;
        EmitVertex();
    }
    EndPrimitive();

}
