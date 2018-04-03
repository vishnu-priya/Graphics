#include "render_world.h"
#include "flat_shader.h"
#include "object.h"
#include "light.h"
#include "ray.h"
#include "sphere.h"


Render_World::Render_World()
  :background_shader(0),ambient_intensity(0),enable_shadows(true),
    recursion_depth_limit(3){}

Render_World::~Render_World()
{
  delete background_shader;
  for(size_t i=0;i<objects.size();i++) delete objects[i];
  for(size_t i=0;i<lights.size();i++) delete lights[i];
}

// Find the closest object of intersection and return the object that was
// intersected.  Record the Hit structure in hit.  If no intersection occurred,
// return NULL.  Note that in the case of a Boolean, the object returned will be
// the Boolean, but the object stored in hit will be the underlying primitive.
// Any intersection with t<=small_t should be ignored.
Object* Render_World::Closest_Intersection(const Ray& ray, Hit& hit)
{
  // TODO
  std::vector<Hit> hits;
  int count = 0;
  hit.object = NULL;
  for(int i =0; i< objects.size(); i++){
    if(objects[i]->Intersection(ray, hits)){
      if(count == 0){
        if(hits.back().ray_exiting){
          hit = hits[hits.size()-1];
        }else{
          hit = hits.back();
        }
        count = 1;
      }else{
        if(hits.back().ray_exiting){
          if(hit.t > hits[hits.size()-1].t){
            hit = hits[hits.size()-1];
          }
        }else{
          if(hit.t > hits.back().t){
            hit = hits.back();
          }
        }
      }
    }
  }
  return 0;
}

// set up the initial view ray and call
void Render_World::Render_Pixel(const ivec2& pixel_index)
{
  Ray ray; // TODO: set up the initial view ray here
  ray.endpoint = camera.position;
  vec3 z = camera.World_Position(pixel_index) - ray.endpoint;
  ray.direction = z/z.magnitude();
  vec3 color=Cast_Ray(ray,1);
  camera.Set_Pixel(pixel_index,Pixel_Color(color));
}

void Render_World::Render()
{
  for(int j=0;j<camera.number_pixels[1];j++)
  for(int i=0;i<camera.number_pixels[0];i++)
  Render_Pixel(ivec2(i,j));
}

// cast ray and return the color of the closest intersected surface point,
// or the background color if there is no object intersection
vec3 Render_World::Cast_Ray(const Ray& ray,int recursion_depth)
{
  // TODO
  vec3 color;
  Hit h;
  Closest_Intersection(ray, h);
  // determine the color here
  if(h.object != NULL){
    vec3 intersection = ray.Point(h.t);
    vec3 normal = h.object->Normal(intersection);
    color = h.object->material_shader->Shade_Surface(ray, intersection, normal, recursion_depth);
  }
  return color;
}
