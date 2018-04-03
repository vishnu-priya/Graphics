#include "light.h"
#include "phong_shader.h"
#include "ray.h"
#include "render_world.h"
#include "object.h"
#include "math.h"


vec3 Phong_Shader::
Shade_Surface(const Ray& ray,const vec3& intersection_point,
    const vec3& same_side_normal,int recursion_depth) const
{
    vec3 color;
    // TODO: determine the color
    vec3 lposition = world.lights[0]->position - intersection_point;
    lposition = lposition/lposition.magnitude();
    vec3 reflectedLight = 2*dot(same_side_normal, lposition)*(same_side_normal - lposition);
    reflectedLight = reflectedLight/reflectedLight.magnitude();
    double diffuseFactor = dot(same_side_normal, lposition);
    if(diffuseFactor < 0){
      diffuseFactor = 0;
    }
    double specularFactor = dot(ray.direction, reflectedLight);
    if(specularFactor < 0){
      specularFactor = 0;
    }
    color = color_ambient + color_diffuse*diffuseFactor + color_specular*pow(specularFactor, specular_power) ;
    return color;
}
