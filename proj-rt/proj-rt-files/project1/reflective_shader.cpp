#include "reflective_shader.h"
#include "ray.h"
#include "render_world.h"
#include "light.h"

vec3 Reflective_Shader::
Shade_Surface(const Ray& ray,const vec3& intersection_point,
    const vec3& same_side_normal,int recursion_depth) const
{
    vec3 color;
    // TODO: determine the color
    for(int i = 0; i< world.lights.size(); i++){
      vec3 lposition = world.lights[i]->position - intersection_point;
      vec3 reflectedLight = 2*dot(same_side_normal, lposition)*(same_side_normal) - lposition;
      reflectedLight = reflectedLight/reflectedLight.magnitude();
      Ray reflectedRay = Ray(intersection_point, reflectedLight);
      vec3 color_surface = shader->Shade_Surface(ray, intersection_point, same_side_normal, recursion_depth);
      color += color_surface + reflectivity*(world.Cast_Ray(reflectedRay, recursion_depth+1) - color_surface);

    }
    return color;
}
