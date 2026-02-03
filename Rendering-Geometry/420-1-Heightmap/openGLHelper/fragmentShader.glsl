#version 150

in vec4 col;
out vec4 c;

void main()
{
  // compute the final pixel color
  c = col;
}



/*

// Fragment shader example
// Assuming 'vColor' is the varying grayscale value from the vertex shader
out vec4 color;
void main() {
    color = vec4(vColor, vColor, vColor, 1.0); // Convert grayscale to RGB
}
*/