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
    color = world.ambient_color*color_ambient*world.ambient_intensity;
    for(int i = 0; i< world.lights.size(); i++){
      vec3 lposition = world.lights[i]->position - intersection_point;
      Ray pointRay = Ray(intersection_point, lposition);
      Hit h;
      double distance = lposition.magnitude();
      world.Closest_Intersection(pointRay, h);
      if(h.t < small_t || h.t > distance || !world.enable_shadows){
        lposition = lposition/distance;
        Ray light_ray = Ray(world.lights[i]->position, lposition);
        vec3 lightEmitted = world.lights[i]->Emitted_Light(light_ray)/(distance*distance);
        vec3 reflectedLight = 2*dot(same_side_normal, lposition)*(same_side_normal) - lposition;
        reflectedLight = reflectedLight/reflectedLight.magnitude();
        double diffuseFactor = dot(same_side_normal, lposition);
        if(diffuseFactor < 0){
          diffuseFactor = 0;
        }
        double specularFactor = dot(-ray.direction, reflectedLight);
        if(specularFactor < 0){
          specularFactor = 0;
        }
        color += lightEmitted * color_diffuse*diffuseFactor + lightEmitted*color_specular*pow(specularFactor, specular_power) ;
    }else{
        color += world.background_shader->Shade_Surface(ray, intersection_point, same_side_normal, recursion_depth);
      }
    }
    return color;
}
