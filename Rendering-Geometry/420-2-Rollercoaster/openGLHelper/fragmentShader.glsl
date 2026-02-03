#version 150

in vec3 viewPosition;
in vec3 viewNormal;
out vec4 c; // output color
uniform vec4 La = vec4(1.0, 1.0, 1.0, 1.0); // light ambient
uniform vec4 Ld = vec4(1.0, 1.0, 1.0, 1.0 ); // light diffuse
uniform vec4 Ls = vec4(1.0, 1.0, 1.0, 1.0 ); // light specular
uniform vec3 viewLightDirection;
uniform vec4 ka = vec4(0.1, 0.1, 0.1, 1.0); // mesh ambient
uniform vec4 kd = vec4(0.1, 0.1, 0.1, 1.0); // mesh diffuse
uniform vec4 ks = vec4(0.2, 0.2, 0.2, 1.0); // mesh specular
uniform float alpha = 1.0; // shininess

void main()
{
    // camera is at (0,0,0) after the modelview transformation
    vec3 eyedir = normalize(vec3(0, 0, 0) - viewPosition);
    // reflected light direction
    vec3 reflectDir = -reflect(viewLightDirection, viewNormal);
    // Phong lighting
    float d = max(dot(viewLightDirection, viewNormal), 0.0f);
    float s = max(dot(reflectDir, eyedir), 0.0f);
    // compute the final color
    c = ka * La + d * kd * Ld + pow(s, alpha) * ks * Ls;
    //c = ka * La + pow(s, alpha) * ks * Ls;
    //c = ka * La ;
    //c = d * kd * Ld;
    //c = pow(s, alpha) * ks * Ls;
}


/*

// Fragment shader example
// Assuming 'vColor' is the varying grayscale value from the vertex shader
out vec4 color;
void main() {
    color = vec4(vColor, vColor, vColor, 1.0); // Convert grayscale to RGB
}
*/