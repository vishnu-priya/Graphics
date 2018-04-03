#include "boolean.h"


// Determine if the ray intersects with the boolean of A and B.
bool Boolean::Intersection(const Ray& ray, std::vector<Hit>& hits) const
{
    // TODO
    std::vector<Hit> hit_calcA, hit_calcB;
    A->Intersection(ray, hit_calcA);
    B->Intersection(ray, hit_calcB);
    A->material_shader = this->material_shader;
    B->material_shader = this->material_shader;
    struct t{
      static void addIntoHits(const Object *o, double t, bool b, std::vector<Hit>& hits){
        Hit h;
        h.object = o;
        h.t = t;
        h.ray_exiting = b;
        hits.push_back(h);
      }
    };
    bool track = 0;
    int i=0, j=0;
    if(hit_calcA.size()||hit_calcB.size()){
      switch (type) {
        case type_union:
          while(i<hit_calcA.size()&&j<hit_calcB.size()){
            if(hit_calcA[i].t < hit_calcB[j].t){
              if(hit_calcB[j].t < hit_calcA[i+1].t){
                t::addIntoHits(hit_calcA[i].object,hit_calcA[i].t,track, hits);
                track = !track;
                if(hit_calcA[i+1].t < hit_calcB[j+1].t){
                  t::addIntoHits(hit_calcB[j+1].object, hit_calcB[j+1].t, track, hits);
                  track = !track;
                }else{
                  t::addIntoHits(hit_calcA[i+1].object, hit_calcA[i+1].t, track, hits);
                  track = !track;
                }
              }else if(hit_calcA[i+1].t < hit_calcB[j].t){
                t::addIntoHits(hit_calcA[i].object,hit_calcA[i].t,0, hits);
                t::addIntoHits(hit_calcA[i+1].object, hit_calcA[i+1].t, 1, hits);
                t::addIntoHits(hit_calcB[j].object, hit_calcB[j].t, 0, hits);
                t::addIntoHits(hit_calcB[j+1].object, hit_calcB[j+1].t, 1, hits);
              }
            }else{
              if(hit_calcA[i].t < hit_calcB[j+1].t){
                t::addIntoHits(hit_calcB[j].object,hit_calcB[j].t,track, hits);
                track = !track;
                if(hit_calcB[i+1].t < hit_calcA[j+1].t){
                  t::addIntoHits(hit_calcA[i+1].object, hit_calcA[i+1].t, track, hits);
                  track = !track;
                }else{
                  t::addIntoHits(hit_calcB[j+1].object, hit_calcB[j+1].t, track, hits);
                  track = !track;
                }
              }else if(hit_calcB[j+1].t < hit_calcA[i].t){
                t::addIntoHits(hit_calcA[i].object,hit_calcA[i].t,0, hits);
                t::addIntoHits(hit_calcA[i+1].object, hit_calcA[i+1].t, 1, hits);
                t::addIntoHits(hit_calcB[j].object, hit_calcB[j].t, 0, hits);
                t::addIntoHits(hit_calcB[j+1].object, hit_calcB[j+1].t, 1, hits);
              }
            }
            i += 2; j+= 2;
          }
          if(i<hit_calcA.size()){
            for(;i<hit_calcA.size(); i++){
              t::addIntoHits(hit_calcA[i].object,hit_calcA[i].t,hit_calcA[i].ray_exiting, hits);
            }
          }
          if(j<hit_calcB.size()){
            for(;j<hit_calcB.size(); j++){
              t::addIntoHits(hit_calcB[j].object,hit_calcB[j].t,hit_calcB[j].ray_exiting, hits);
            }
          }
          break;
        case type_intersection:
          while(i<hit_calcA.size()&&j<hit_calcB.size()){
            if(hit_calcA[i].t < hit_calcB[j].t){
              if(hit_calcB[j].t < hit_calcA[i+1].t){
                t::addIntoHits(hit_calcB[j].object,hit_calcB[j].t,0, hits);
                if(hit_calcA[i+1].t < hit_calcB[j+1].t){
                  t::addIntoHits(hit_calcA[i+1].object, hit_calcA[i+1].t, 1, hits);
                }else{
                  t::addIntoHits(hit_calcB[j+1].object, hit_calcB[j+1].t, 1, hits);
                }
              }else if(hit_calcB[j].t > hit_calcA[i+1].t){
                return false;
              }
            }else{
              if(hit_calcA[i].t < hit_calcB[j+1].t){
                t::addIntoHits(hit_calcA[i].object,hit_calcA[i].t,0, hits);
                if(hit_calcB[j+1].t < hit_calcA[i+1].t){
                  t::addIntoHits(hit_calcB[j+1].object, hit_calcB[j+1].t, 1, hits);
                }else{
                  t::addIntoHits(hit_calcA[i+1].object, hit_calcA[i+1].t, 1, hits);
                }
              }else if(hit_calcA[j].t > hit_calcB[i+1].t){
                return false;
              }
            }
            i+=2; j+=2;
          }
          if(i<hit_calcA.size()||j<hit_calcB.size()){
            return false;
          }
          break;
        case type_difference:
          while(i<hit_calcA.size()&&j<hit_calcB.size()){
            for(int a = 0; a< hit_calcA.size(); a++){
            }
            for(int a = 0; a< hit_calcB.size(); a++){
            }
            if(hit_calcA[i].t < hit_calcB[j].t){
              t::addIntoHits(hit_calcA[i].object,hit_calcA[i].t,0, hits);
              if(hit_calcB[j].t < hit_calcA[i+1].t){
                if(hit_calcA[i+1].t < hit_calcB[j+1].t){
                  t::addIntoHits(hit_calcB[j].object, hit_calcB[j].t, 1, hits);
                }else{
                  t::addIntoHits(hit_calcB[j].object, hit_calcB[j].t, 1, hits);
                  t::addIntoHits(hit_calcA[i+1].object, hit_calcA[i+1].t, 0, hits);
                  t::addIntoHits(hit_calcB[j+1].object, hit_calcB[j+1].t, 1, hits);
                }
              }else if(hit_calcB[j].t > hit_calcA[i+1].t){
                t::addIntoHits(hit_calcA[i].object,hit_calcA[i].t,0, hits);
                t::addIntoHits(hit_calcA[i+1].object, hit_calcA[i+1].t, 0, hits);
              }
            }else{
              if(hit_calcA[i].t < hit_calcB[j+1].t){
                if(hit_calcA[i+1].t < hit_calcB[j+1].t){
                  return false;
                }
                t::addIntoHits(hit_calcB[j+1].object,hit_calcB[j+1].t,0, hits);
                t::addIntoHits(hit_calcA[i+1].object, hit_calcA[i+1].t, 1, hits);
              }else{
                t::addIntoHits(hit_calcA[i].object, hit_calcA[i].t, 1, hits);
                t::addIntoHits(hit_calcA[i+1].object, hit_calcA[i+1].t, 1, hits);
              }
            }
            i+=2; j+=2;
          }
          if(i<hit_calcA.size()){
            for(;i<hit_calcA.size(); i++){
              t::addIntoHits(hit_calcA[i].object,hit_calcA[i].t,hit_calcA[i].ray_exiting, hits);
            }
          }
          if(j<hit_calcB.size()){
            return false;
          }
          break;
      }
      return true;
    }
    return false;
}

// This should never be called.
vec3 Boolean::Normal(const vec3& point) const
{
    assert(false);
    return vec3();
}
